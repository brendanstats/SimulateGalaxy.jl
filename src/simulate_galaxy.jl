"""
Function to simulate a galaxy depending on type of parameters specified given the sample size.
"""
function simulate_galaxy{G <: AbstractFloat, T <: Integer}(p::NFWParameters{G}, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false)
    function nfwopt(x::Array{G, 1})
        return -(x[1] * x[1] * x[3] * profile_density(x[1], x[2], x[3], p, αt = αt))
    end

    minProb = Optim.optimize(nfwopt, [0.1, 1.0, 1.0], method = Optim.NelderMead())
    fmax = max(0.0,  -1.1 * minProb.minimum)
    
    vmax0 = escape_velocity(x0, p)
    sampledr = Array{Float64}(samplesize)
    sampledvr = Array{Float64}(samplesize)
    sampledvt = Array{Float64}(samplesize)
    if rate
        tested = 0
    end
    accepted = 0
    while accepted < samplesize
        x = (p.rlim  / p.rs) * rand()
        vr = vmax0 * rand()
        vt = vmax0 * rand()
        u = rand()
        sampleDensity = x * x * vt * profile_density(x, vr, vt, p, αt = αt)
        if rate
            tested += 1
        end
        if sampleDensity >= (u * fmax)
            accepted += 1
            sampledr[accepted] = x * p.rs
            sampledvr[accepted] =  vr
            sampledvt[accepted] = vt
        end
    end
    if rate
        return SphericalGalaxy(sampledr, sampledvr, sampledvt, samplesize), accepted / tested
    end
    return SphericalGalaxy(sampledr, sampledvr, sampledvt, samplesize)
end

function simulate_galaxy{G <: AbstractFloat, T <: Integer}(p::SFWParameters{G}, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false)
    function sfwopt(x::Array{G, 1})
        return -(x[1] * x[1] * x[3] * profile_density(x[1], x[2], x[3], p, αt = αt))
    end

    minProb = Optim.optimize(sfwopt, [0.1, 1.0, 1.0], method = Optim.NelderMead())
    fmax = max(0.0,  -1.1 * minProb.minimum)
    
    vmax0 = escape_velocity(x0, p)
    sampledr = Array{Float64}(samplesize)
    sampledvr = Array{Float64}(samplesize)
    sampledvt = Array{Float64}(samplesize)
    if rate
        tested = 0
    end
    accepted = 0
    while accepted < samplesize
        x = (p.rlim  / p.rs) * rand()
        vr = vmax0 * rand()
        vt = vmax0 * rand()
        u = rand()
        sampleDensity = x * x * vt * profile_density(x, vr, vt, p, αt = αt)
        if rate
            tested += 1
        end
        if sampleDensity >= (u * fmax)
            accepted += 1
            sampledr[accepted] = x * p.rs
            sampledvr[accepted] =  vr
            sampledvt[accepted] = vt
        end
    end
    if rate
        return SphericalGalaxy(sampledr, sampledvr, sampledvt, samplesize), accepted / tested
    end
    return SphericalGalaxy(sampledr, sampledvr, sampledvt, samplesize)
end


function simulate_galaxy{G <: AbstractFloat, T <: Integer}(p::NFWParameters{G}, mμ::G, mσ::G, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false)
    metal = randn(samplesize) .* mσ .+ mμ
    if rate
        sg, acptrt = simulate_galaxy(p, samplesize, x0 = x0, αt = αt, rate = rate)
        return MetallicSphericalGalaxy(sg, metal), acptrt
    else
        return MetallicSphericalGalaxy(simulate_galaxy(p, samplesize, x0 = x0, αt = αt, rate = rate), metal)
    end
end

function simulate_galaxy{G <: AbstractFloat, T <: Integer}(p::SFWParameters{G}, mμ::G, mσ::G, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false)
    metal = randn(samplesize) .* mσ .+ mμ
    if rate
        sg, acptrt = simulate_galaxy(p, samplesize, x0 = x0, αt = αt, rate = rate)
        return MetallicSphericalGalaxy(sg, metal), acptrt
    else
        return MetallicSphericalGalaxy(simulate_galaxy(p, samplesize, x0 = x0, αt = αt, rate = rate), metal)
    end
end
