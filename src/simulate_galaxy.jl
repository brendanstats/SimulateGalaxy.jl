"""
Function to simulate a galaxy depending on type of parameters specified given the sample size.
"""
function simulate_galaxy(p::NFWParameters{G}, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false) where {G <: AbstractFloat, T <: Integer}
    function nfwopt(x::Array{G, 1})
        return -(x[1] * x[1] * x[3] * profile_density(x[1], x[2], x[3], p, αt = αt))
    end
    #function nfwopt(x::Array{G, 1})
    #    return -(2.0 * log(x[1]) +  log(x[3]) + log_profile_density(x[1], x[2], x[3], p, αt = αt))
    #end

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
        #u = log(rand())
        sampleDensity = x * x * vt * profile_density(x, vr, vt, p, αt = αt)
        #sampleDensity = 2.0 * log(x) + log(vt) + log_profile_density(x, vr, vt, p, αt = αt)
        if rate
            tested += 1
        end
        if sampleDensity >= (u * fmax)
            #if sampleDensity >= (u + fmax)
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

function simulate_galaxy(p::SFWParameters{G}, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false) where {G <: AbstractFloat, T <: Integer}
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


function simulate_galaxy(p::NFWParameters{G}, mμ::G, mσ::G, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false) where {G <: AbstractFloat, T <: Integer}
    metal = randn(samplesize) .* mσ .+ mμ
    if rate
        sg, acptrt = simulate_galaxy(p, samplesize, x0 = x0, αt = αt, rate = rate)
        return MetallicSphericalGalaxy(sg, metal), acptrt
    else
        return MetallicSphericalGalaxy(simulate_galaxy(p, samplesize, x0 = x0, αt = αt, rate = rate), metal)
    end
end

function simulate_galaxy(p::SFWParameters{G}, mμ::G, mσ::G, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false) where {G <: AbstractFloat, T <: Integer}
    metal = randn(samplesize) .* mσ .+ mμ
    if rate
        sg, acptrt = simulate_galaxy(p, samplesize, x0 = x0, αt = αt, rate = rate)
        return MetallicSphericalGalaxy(sg, metal), acptrt
    else
        return MetallicSphericalGalaxy(simulate_galaxy(p, samplesize, x0 = x0, αt = αt, rate = rate), metal)
    end
end


"""
Use a grid rejection sampler instead of a standard rejection sampler
"""
function simulate_galaxy_grid(p::SFWParameters{G}, samplesize::T; x0::G = 10.0^-8, αt::G = 1.0, rate::Bool = false, ngrid::Integer = 3) where {G <: AbstractFloat, T <: Integer}
    
    vmax0 = escape_velocity(x0, p)
    
    function sfwopt(x::Array{G, 1})
        return -(x[1] * x[1] * x[3] * profile_density(x[1], x[2], x[3], p, αt = αt))
    end

    ##grid bounds
    xgrid = linspace(0.0, p.rlim  / p.rs, ngrid + 1)
    vrgrid = linspace(0.0, vmax0, ngrid + 1)
    vtgrid = linspace(0.0, vmax0, ngrid + 1)

    ##midpoints
    xmids = xgrid[1:ngrid] + 0.5 * step(xgrid)
    vrmids = vrgrid[1:ngrid] + 0.5 * step(vrgrid)
    vtmids = vtgrid[1:ngrid] + 0.5 * step(vtgrid)

    ##optimization variables
    startpoints = Array{Float64}(ngrid^3, 3)
    lowers = Array{Float64}(ngrid^3, 3)
    uppers = Array{Float64}(ngrid^3, 3)
    fmaxes = Array{Float64}(ngrid^3)
    ii = 1
    for xx in 1:ngrid, rr in 1:ngrid, tt in 1:ngrid
        startpoints[ii, :] = [xmids[xx], vrmids[rr], vtmids[tt]]
        lowers[ii, :] = [xgrid[xx], vrgrid[rr], vtgrid[tt]]
        uppers[ii, :] = [xgrid[xx + 1], vrgrid[rr + 1], vtgrid[tt + 1]]
        minProb = Optim.optimize(sfwopt, startpoints[ii, :], lowers[ii, :], uppers[ii, :], Fminbox{Optim.NelderMead}())
        fmaxes[ii] = max(0.0,  -1.1 * minProb.minimum)
        ii += 1
    end
    w = StatsBase.Weights(fmaxes)  
    #http://julianlsolvers.github.io/Optim.jl/latest/user/minimization/#box-minimization
        
    sampledr = Array{Float64}(samplesize)
    sampledvr = Array{Float64}(samplesize)
    sampledvt = Array{Float64}(samplesize)
    if rate
        tested = 0
    end
    accepted = 0
    while accepted < samplesize
        idx = StatsBase.sample(w)
        x = lowers[ii, 1] + step(xgrid) * rand()
        vr = lowers[ii, 2] + step(vrgrid) * rand()
        vt = lowers[ii, 3] + step(vtgrid) * rand()
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
