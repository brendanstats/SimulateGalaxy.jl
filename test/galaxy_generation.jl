#=
Sampling functions to generate random galaxy consistent with provided NFW parameters
=#

#=
Questions from while coding
1) Why more complicated PhiR (genphi()) function in previous functions
2) Why x^2*vt in fprob()
3) Why introduce dummy normalization? in f()
=#

using Optim, GSL

#=
Define gravitationla potential function,
results are in units of PhiS
=#
function GravPotential(x::Real)
    return 1 - (log(1 + x) / x)
end

#=
Calculate escape velocity for a radius scaled by the
characteristic radius, i.e. x = r / rs
=#
function EscapeVelocity(x, rlim, rs, PhiS)
    xlim = rlim / rs
    return sqrt(2 * PhiS * (GravPotential(xlim) - GravPotential(x)))
end

#=
Calculate the  Navarro-Frenk-White density based on x, the radius
scaled by the characteristic radius, vr the radial component of velocity,
vt the tangental component of velocity, and the NFW parameters.  The density
is not normalized, the normalization constant is set to 1000 for conveniance.
=#
function NFWDensity{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1})
    if x < 0
        return 0
    end
    #=
    a = parameters[1]
    d = parameters[2]
    e = parameters[3]
    Ec = parameters[4]
    rlim = parameters[5]
    b = parameters[6]
    q = parameters[7]
    Jb = parameters[8]
    vmax = parameters[9]
    rmax = parameters[10]
    =#

    a, d, e, Ec, rlim, b, q, Jb, vmax, rmax = parameters

    rs = rmax / 2.16
    PhiS = (vmax / 0.465)^2
    PhiR = GravPotential(x) * PhiS

    J = abs(x * rs * vt) #vt = v * sin(theta)
    E = (vt * vt + vr * vr) / 2 + PhiR

    Ec = Ec * PhiS
    xlim = rlim / rs
    Philim = PhiS * GravPotential(xlim)
    Jb = Jb * rs * sqrt(PhiS)

    if b <= 0
        gJ = 1 / (1 + (J / Jb)^(-b))
    else
        gJ = 1 + (J / Jb)^b
    end

    N = 1000
    if E < Philim && E >= 0
        hE = N * E^a * (E^q + Ec^q)^(d/q) * (Philim - E)^e
    else
        hE = 0
    end
    return hE * gJ
end

#=
Scale NFWDensity function
=#
function NFWProb{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1})
    return x * x * vt * NFWDensity(x, vr, vt, parameters)
end

#=
Randomly a generate with samplesize = # stars and NFW parameters parameters.
Basic rejection sampler with a factor of 1.1 is used.  Results are scaled radius (x),
radial velocity (vr), and tangential velocity (vt).
=#
function GenGalaxy{G<:Real}(parameters::Array{G,1}, samplesize::Integer)
    rlim = parameters[5]
    vmax = parameters[9]
    rmax = parameters[10]
    rs = rmax / 2.16
    PhiS = (vmax / 0.465)^2
    
    function NFWProbOpt(x)
        return -NFWProb(x[1], x[2], x[3], parameters)
    end

    minProb = optimize(NFWProbOpt, [0.1, 1.0, 1.0], method = :nelder_mead)
    fmax = max(0,  -1.1 * minProb.f_minimum)
    x0 = 10.0^-8
    vmax0 = EscapeVelocity(x0, rlim, rs, PhiS)
    #print("fmax: ", fmax, " vmax: ", vmax0, "\n")

    sampledValues = Array{Array{Float64,1}}(samplesize)
    tested = 0
    accepted = 0
    while accepted < samplesize
        x = (rlim  / rs) * rand()
        vr = vmax0 * rand()
        vt = vmax0 * rand()
        u = rand()
        sampleDensity = NFWProb(x, vr, vt, parameters)
        tested += 1
        if sampleDensity >= (u * fmax)
            accepted += 1
            sampledValues[accepted] = [x*rs , vr, vt]
        end
    end
    #print("acceptance rate:", accepted / tested, "\n")
    return sampledValues
end

#=
Transform the results of GenGalaxy to Euclidean coordinates.
=#
function converttoeuclidean{G<:Real}(coordinate::Array{G,1})
    radius = coordinate[1]
    theta = acos(1.0 - 2.0 * rand())
    phi = 2.0 * pi * rand()

    velocity = sqrt(coordinate[2]^2 + coordinate[3]^2)
    velSign = sign(rand() - 0.5)
    velPhi = 2.0 * pi * rand()

    x = radius * sin(theta) * cos(phi)
    y = radius * sin(theta) * sin(phi)
    z = radius * cos(theta)

    vz2 = velSign * coordinate[2]
    vx2 = coordinate[3] * cos(velPhi)
    vy2 = coordinate[3] * sin(velPhi)

    vx = cos(theta) * cos(phi) * vx2 - sin(phi) * vy2 + sin(theta) * cos(phi) * vz2
    vy = cos(theta) * sin(phi) * vx2 + cos(phi) * vy2 + sin(theta) * sin(phi) * vz2
    vz = -sin(theta) * vx2 + cos(theta) * vz2
    return x, y, z, vx, vy, vz
end

#=
Simplified posterior functions
=#
function NFWhEDensity{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1})
    a = parameters[1]
    d = parameters[2]
    e = parameters[3]
    Ec = parameters[4]
    rlim = parameters[5]
    b = parameters[6]
    q = parameters[7]
    Jb = parameters[8]
    vmax = parameters[9]
    rmax = parameters[10]

    rs = rmax / 2.16
    PhiS = (vmax / 0.465)^2
    PhiR = GravPotential(x) * PhiS

    J = abs(x * rs * vt) #vt = v * sin(theta)
    E = (vt * vt + vr * vr) / 2 + PhiR

    Ec = Ec * PhiS
    xlim = rlim / rs
    Philim = PhiS * GravPotential(xlim)
    Jb = Jb * rs * sqrt(PhiS)

    N = 1000
    if E < Philim && E >= 0
        hE = N * E^a * (E^q + Ec^q)^(d/q) * (Philim - E)^e
    else
        hE = 0
    end
    return hE
end

function NFWhEProb{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1})
    return x * x * vt * NFWhEDensity(x, vr, vt, parameters)
end


function GenGalaxySimplified{G<:Real}(parameters::Array{G,1}, samplesize::Integer)
    rlim = parameters[5]
    vmax = parameters[9]
    rmax = parameters[10]
    rs = rmax / 2.16
    PhiS = (vmax / 0.465)^2
    
    function NFWhEProbOpt(x)
        return -NFWhEProb(x[1], x[2], x[3], parameters)
    end

    minProb = optimize(NFWhEProbOpt, [0.1, 1.0, 1.0], method = :nelder_mead)
    fmax = max(0,  -1.1 * minProb.f_minimum)

    sampledValues = Array{Array{Float64,1}}(samplesize)
    tested = 0
    accepted = 0
    while accepted < samplesize
        x = (rlim  / rs) * rand()
        x0 = 10.0^-8
        vmax0 = EscapeVelocity(x0, rlim, rs, PhiS)
        vr = vmax0 * rand()
        vt = vmax0 * rand()
        u = rand()
        sampleDensity = NFWhEProb(x, vr, vt, parameters)
        tested += 1
        if sampleDensity >= (u * fmax)
            accepted += 1
            sampledValues[accepted] = [x*rs , vr, vt]
        end
    end
    #print("acceptance rate:", accepted / tested, "\n")
    return sampledValues
end


function convertgalaxy(galaxy::Array{Array{Float64, 1}, 1})
    n = length(galaxy)
    x = Array{Float64}(n)
    y = Array{Float64}(n)
    z = Array{Float64}(n)
    vx = Array{Float64}(n)
    vy = Array{Float64}(n)
    vz = Array{Float64}(n)
    for ii in 1:n
        x[ii], y[ii], z[ii], vx[ii], vy[ii], vz[ii] = converttoeuclidean(galaxy[ii])
    end
    return x, y, z, vx, vy, vz
end

function radiusvelocity(galaxy::Array{Array{Float64, 1}, 1})
    n = length(galaxy)
    r = Array{Float64}(n)
    v = Array{Float64}(n)
    for ii in 1:n
        r[ii] = galaxy[ii][1]
        v[ii] = sqrt(galaxy[ii][2]^2 + galaxy[ii][3]^2)
    end
    return r, v
end

function radiusvelocity(x::Array{Float64}, y::Array{Float64}, z::Array{Float64}, vx::Array{Float64}, vy::Array{Float64}, vz::Array{Float64})
    n = length(x)
    r = Array{Float64}(n)
    v = Array{Float64}(n)
    for ii in 1:n
        r[ii] = sqrt(x[ii]^2 + y[ii]^2 + z[ii]^2)
        v[ii] = sqrt(vx[ii]^2 + vy[ii]^2 + vz[ii]^2)
    end
    return r, v
end

function radius(x::Array{Float64}, y::Array{Float64})
    n = length(x)
    r = Array{Float64}(n)
    for ii in 1:n
        r[ii] = sqrt(x[ii]^2 + y[ii]^2)
    end
    return r
end

function rv2(r::Array{Float64, 1}, v::Array{Float64, 1})
    n = length(r)
    quantity = Array{Float64}(n)
    for ii in 1:n
        quantity[ii] = r[ii] * v[ii]^2
    end
    return quantity
end

#new functions

function genphi0_integral(x::Float64, alpha::Float64, beta::Float64, gamma::Float64, rhos::Float64, rs::Float64)
    PhiS = rhos * rs^2
    function integrand1(xi)
        return (xi ^ -gamma) * (1 + xi ^ alpha) ^ ((gamma - beta) / alpha) * xi * xi
    end
    function integrand2(xi)
        return (xi ^ -gamma) * (1 + xi ^ alpha) ^ ((gamma - beta) / alpha) * xi
    end
    phi1 = quadgk(integrand1, 0.0, x)[1]
    phi2 = quadgk(integrand2, x, Inf)[1]
    
    return PhiS * (1 - (phi1 / x + phi2))
end

function hypergeom_2F1_expanded(a::Float64, b::Float64, c::Float64, z::Float64)
    #if z >= -1.0 && z < 1.0
    #    return hypergeom([a, b], c, z)
    #end
    if z >= -1.0 && z < 0.5
        return hypergeom([a, b], c, z)
    end
    if z >= 0.5 && z < 1.0
        coef1 = gamma(c) * gamma(c - a - b) / (gamma(c - a) * gamma(c - b))
        hypergeom1 = hypergeom_2F1_expanded(a, b, a + b - c + 1, 1 - z)
        coef2 = (1 - z) ^ (c - a - b) * gamma(c) * gamma(a + b - c) / (gamma(a) * gamma(b))
        hypergeom2 = hypergeom_2F1_expanded(c - a, c - b, c - a - b + 1, 1 - z)
        return coef1 * hypergeom1 + coef2 * hypergeom2
    end
    if z >= 1.0
        coef1 = z ^ (-1.0 * a) * gamma(c) * gamma(c - a - b) / (gamma(c - a) * gamma(c - b))
        hypergeom1 = hypergeom_2F1_expanded(a, a - c + 1, a + b - c + 1, 1 - 1 / z)
        coef2 = z ^ (a - c) * (1 - z) ^ (c - a - b) * gamma(c) * gamma(a + b - c) / (gamma(a) * gamma(b))
        hypergeom2 = hypergeom_2F1_expanded(c - a, 1 - a, c - a - b + 1, 1 - 1 / z)
        return coef1 * hypergeom1 + coef2 * hypergeom2
    end
    if z < -1.0
        return (1 - z) ^ (-1.0 * a) * hypergeom_2F1_expanded(a, c - b, c, z / (z - 1))
    end
end


function genphi0(x::Float64, alpha::Float64, beta::Float64, gamma::Float64, rhos::Float64, rs::Float64)
    PhiS = rhos * rs^2
    x0 = 10.0 ^ -13
    p1a = 0
    p1b = hypergeom_2F1_expanded((3.0 - gamma) / alpha, (beta - gamma) / alpha, (3.0 + alpha - gamma) / alpha, -x ^ alpha)
    I1 = (x0 ^ (3.0 - gamma) * p1a - x ^ (3.0 - gamma) * p1b) / (x * (gamma - 3.0))
    p2 = hypergeom_2F1_expanded((beta - 2.0) / alpha, (beta - gamma) / alpha, (alpha + beta - 2.0) / alpha, -x ^ (-alpha))
    I2 = x ^ (2.0 - beta) * p2 / (beta - 2.0)
    return PhiS * (1 - (I1 + I2))
end


function genphi(x::Float64, alpha::Float64, beta::Float64, gamma::Float64, rhos::Float64, rs::Float64)
    x0 = 10.0 ^ -8
    return genphi0(x, alpha, beta, gamma, rhos, rs) -  genphi0(x0, alpha, beta, gamma, rhos, rs)
end

function genphi(x::Float64, alpha::Float64, beta::Float64, gamma::Float64, rhos::Float64, rs::Float64, Phi0::Float64)
    return genphi0(x, alpha, beta, gamma, rhos, rs) -  Phi0
end


#new vesc function line 385
function escape_velocity(x, rlim, rhos, rs, alpha, beta, gamma)
    xlim = rlim / rs
    PhiR = genphi(x, alpha, beta, gamma, rhos, rs)
    Philim = genphi(xlim, alpha, beta, gamma, rhos, rs)
    return sqrt(2 * (Philim - PhiR))
end

function escape_velocity(x, rlim, rhos, rs, alpha, beta, gamma, Philim)
    PhiR = genphi(x, alpha, beta, gamma, rhos, rs)
    return sqrt(2 * (Philim - PhiR))
end

function escape_velocity(x, rlim, rhos, rs, alpha, beta, gamma, Philim, Phi0)
    PhiR = genphi(x, alpha, beta, gamma, rhos, rs, Phi0)
    return sqrt(2 * (Philim - PhiR))
end

#new f function
function SFW_density{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1})
    if x < 0
        return 0
    end
    a, d, e, Ec, rlim, b, q, Jb, rhos, rs, alpha, beta, gamma = parameters

    PhiS = rhos * rs ^ 2
    PhiR = genphi(x, alpha, beta, gamma, rhos, rs)

    J = abs(x * rs * vt) #vt = v * sin(theta)
    E = (vt * vt + vr * vr) / 2.0 + PhiR

    Ec = Ec * PhiS
    xlim = rlim / rs
    Philim = genphi(xlim, alpha, beta, gamma, rhos, rs)
    Jb = Jb * rs * sqrt(PhiS)

    if b <= 0
        gJ = 1.0 / (1.0 + (J / Jb)^(-b))
    else
        gJ = 1.0 + (J / Jb)^b
    end

    N = 1000.0
    if E < Philim && E >= 0
        hE = N * E^a * (E^q + Ec^q)^(d/q) * (Philim - E)^e
    else
        hE = 0
    end
    return hE * gJ
end

function SFW_density{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1}, Philim::G)
    if x < 0
        return 0
    end
    a, d, e, Ec, rlim, b, q, Jb, rhos, rs, alpha, beta, gamma = parameters

    PhiS = rhos * rs ^ 2
    PhiR = genphi(x, alpha, beta, gamma, rhos, rs)

    J = abs(x * rs * vt) #vt = v * sin(theta)
    E = (vt * vt + vr * vr) / 2.0 + PhiR

    Ec = Ec * PhiS
    Jb = Jb * rs * sqrt(PhiS)

    if b <= 0
        gJ = 1.0 / (1.0 + (J / Jb)^(-b))
    else
        gJ = 1.0 + (J / Jb)^b
    end

    N = 1000.0
    if E < Philim && E >= 0
        hE = N * E^a * (E^q + Ec^q)^(d/q) * (Philim - E)^e
    else
        hE = 0
    end
    return hE * gJ
end

function SFW_density{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1}, Philim::G, Phi0::G)
    if x < 0
        return 0
    end
    a, d, e, Ec, rlim, b, q, Jb, rhos, rs, alpha, beta, gamma = parameters

    PhiS = rhos * rs ^ 2
    PhiR = genphi(x, alpha, beta, gamma, rhos, rs, Phi0)

    J = abs(x * rs * vt) #vt = v * sin(theta)
    E = (vt * vt + vr * vr) / 2.0 + PhiR

    Ec = Ec * PhiS
    Jb = Jb * rs * sqrt(PhiS)

    if b <= 0
        gJ = 1.0 / (1.0 + (J / Jb)^(-b))
    else
        gJ = 1.0 + (J / Jb)^b
    end

    N = 1000.0
    if E < Philim && E >= 0
        hE = N * E^a * (E^q + Ec^q)^(d/q) * (Philim - E)^e
    else
        hE = 0
    end
    return hE * gJ
end


#new fprob
function SFW_prob{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1})
    return x * x * vt * SFW_density(x, vr, vt, parameters)
end

function SFW_prob{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1}, Philim)
    return x * x * vt * SFW_density(x, vr, vt, parameters, Philim)
end

function SFW_prob{G<:Real}(x::G, vr::G, vt::G, parameters::Array{G,1}, Philim, Phi0)
    return x * x * vt * SFW_density(x, vr, vt, parameters, Philim, Phi0)
end


#new sample function
function general_gen_galaxy{G<:Real}(parameters::Array{G,1}, samplesize::Integer)
    a, d, e, Ec, rlim, b, q, Jb, rhos, rs, alpha, beta, gamma = parameters
    xlim = rlim / rs
    x0 = 10.0 ^ -8
    Phi0 = genphi0(x0, alpha, beta, gamma, rhos, rs)
    Philim = genphi(xlim, alpha, beta, gamma, rhos, rs, Phi0)
    function SFW_prob_opt(x)
        return -SFW_prob(x[1], x[2], x[3], parameters, Philim, Phi0)
    end

    minProb = optimize(SFW_prob_opt, [0.1, 1.0, 1.0], method = :nelder_mead)
    #minProb = optimize(SFW_prob_opt, [0.1, 1.0, 1.0], method = Optim.NelderMead())
    fmax = max(0,  -1.1 * minProb.f_minimum)
    x0 = 10.0^-8
    vmax0 = escape_velocity(x0, rlim, rhos, rs, alpha, beta, gamma, Philim, Phi0)
    #print("fmax: ", fmax, " vmax: ", vmax0, "\n")

    sampledValues = Array{Array{Float64,1}}(samplesize)
    tested = 0
    accepted = 0
    while accepted < samplesize
        x = (rlim  / rs) * rand()
        vr = vmax0 * rand()
        vt = vmax0 * rand()
        u = rand()
        sampleDensity = SFW_prob(x, vr, vt, parameters, Philim, Phi0)
        tested += 1
        if sampleDensity >= (u * fmax)
            accepted += 1
            sampledValues[accepted] = [x * rs , vr, vt]
        end
    end
    #print("acceptance rate:", accepted / tested, "\n")
    return sampledValues
end

#=
param=[2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 4229.2, .694444, 1.01, 3.01, 1.01]
@time general_gen_galaxy(param, 1000); #7.5 seconds before Philim
@time general_gen_galaxy(param, 1000); #3.6 seconds before Phi0
@time general_gen_galaxy(param, 1000); #2.1 seconds with Philim and Phi0

#hypergeom([1.0, 2.0], 2.0, -0.98) #julia syntax
#ss.hyp2f1(1.0, 2.0, 2.0, -1.0) #python syntax
hypergeom()

@time genphi0(x, alpha, beta, gamma, rhos, rs)
@time genphi0_2(1.1, alpha, beta, gamma, rhos, rs)

@everywhere function f(n::Int64)
    x, vt, vr = general_gen_galaxy(param, 100)
    return mean(x)
end
=#
