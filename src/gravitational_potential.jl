"""
More general version of hypergeometric2F1 function than is provided by GSL library, based on 'Computation of Hypergeometric Functions' (Pearson 2009)
"""
function hypergeom_2F1_expanded{G <: AbstractFloat}(a::G, b::G, c::G, z::G)
    if z < -1.0
        return (1.0 - z) ^ (-1.0 * a) * hypergeom_2F1_expanded(a, c - b, c, z / (z - 1.0))
    elseif z < 0.5
        return hypergeom([a, b], c, z)
    elseif z < 1.0
        coef1 = gamma(c) * gamma(c - a - b) / (gamma(c - a) * gamma(c - b))
        hypergeom1 = hypergeom_2F1_expanded(a, b, a + b - c + 1.0, 1.0 - z)
        coef2 = (1.0 - z) ^ (c - a - b) * gamma(c) * gamma(a + b - c) / (gamma(a) * gamma(b))
        hypergeom2 = hypergeom_2F1_expanded(c - a, c - b, c - a - b + 1.0, 1.0 - z)
        return coef1 * hypergeom1 + coef2 * hypergeom2
    else #if z >= 1.0
        coef1 = z ^ (-1.0 * a) * gamma(c) * gamma(c - a - b) / (gamma(c - a) * gamma(c - b))
        hypergeom1 = hypergeom_2F1_expanded(a, a - c + 1.0, a + b - c + 1.0, 1.0 - 1.0 / z)
        coef2 = z ^ (a - c) * (1.0 - z) ^ (c - a - b) * gamma(c) * gamma(a + b - c) / (gamma(a) * gamma(b))
        hypergeom2 = hypergeom_2F1_expanded(c - a, 1.0 - a, c - a - b + 1.0, 1.0 - 1.0 / z)
        return coef1 * hypergeom1 + coef2 * hypergeom2
    end
end


"""
Dark matter Density profile given in units of ρs
"""
function density_profile{G <: AbstractFloat}(x::G)
    return (x * (1.0 + x) ^ 2) ^ -1
end

function density_profile{G <: AbstractFloat, P <: SFWParameters}(x::G, p::P)
    return (x ^ -p.γ) * (1.0 + x ^ p.α) ^ ((p.γ - p.β) / p.α)
end


"""
Define gravitationla potential function, results are in units of Φs (i.e. are often multiplied by it)
"""
function gravitational_potential{G <: AbstractFloat}(x::G)
    return 1.0 - (log(1.0 + x) / x)
end

function gravitational_potential{G <: AbstractFloat}(x::G, p::NFWParameters{G})
    return p.Φs * (1.0 - (log(1.0 + x) / x))
end

"""
Slower version of function
"""
function gravitational_potential_integral{G <: AbstractFloat, P <: SFWParameters}(x::G, p::P)
    function integrand1(xi::G)
        return  density_profile(xi, p) * xi * xi
    end
    function integrand2(xi::G)
        return density_profile(xi, p) * xi
    end
    Φ1 = quadgk(integrand1, 0.0, x)[1]
    Φ2 = quadgk(integrand2, x, Inf)[1]
    
    return p.Φs * (1.0 - (Φ1 / x + Φ2))
end

function gravitational_potential{G <: AbstractFloat, P <: SFWParameters}(x::G, p::P)
    if ((3.0 - p.β) / p.α % 1.0) == 0.0 || ((p.γ - p.β) / p.α % 1.0) == 0.0
        return gravitational_potential_integral(x, p)
    end
    p1b = hypergeom_2F1_expanded((3.0 - p.γ) / p.α, (p.β - p.γ) / p.α, (3.0 + p.α - p.γ) / p.α, -x ^ p.α)
    I1 = -(x ^ (3.0 - p.γ) * p1b) / (x * (p.γ - 3.0))
    p2 = hypergeom_2F1_expanded((p.β - 2.0) / p.α, (p.β - p.γ) / p.α, (p.α + p.β - 2.0) / alpha, -x ^ (-p.α))
    I2 = x ^ (2.0 - p.β) * p2 / (p.β - 2.0)
    return p.Φs * (1.0 - (I1 + I2))
end

"""
Gravitational potential for SFW profile where center of galaxy is forced to have a potential of 0.
"""
function gravitational_potential0{G <: AbstractFloat, P <: SFWParameters}(x::G, p::P)
    return gravitational_potential(x, p) -  p.Φ0
end
