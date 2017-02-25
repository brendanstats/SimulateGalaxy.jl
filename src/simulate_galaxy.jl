"""
Function to simulate a galaxy depending on type of parameters specified given the sample size.
"""
function simulate_galaxy{G <: AbstractFloat, T <: Integer}(p::NFWParameters{G}, samplesize::T; x0::G = 10.0^-8, rate::Bool = false)
    function nfwopt(x::Array{G, 1})
        return -(x[1] * x[1] * x[3] * profile_density(x[1], x[2], x[3], p))
    end

    minProb = Optim.optimize(nfwopt, [0.1, 1.0, 1.0], method = :nelder_mead)
    fmax = max(0.0,  -1.1 * minProb.f_minimum)
    
    vmax0 = escape_velocity(x0, p)
    sampledValues = Array{Float64}(samplesize, 3)
    if rate
        tested = 0
    end
    accepted = 0
    while accepted < samplesize
        x = (p.rlim  / p.rs) * rand()
        vr = vmax0 * rand()
        vt = vmax0 * rand()
        u = rand()
        sampleDensity = x * x * vt * profile_density(x, vr, vt, p)
        if rate
            tested += 1
        end
        if sampleDensity >= (u * fmax)
            accepted += 1
            sampledValues[accepted, :] = [x * p.rs , vr, vt]
        end
    end
    if rate
        return sampledValues, accepted / tested
    end
    return sampledVaues
end

function simulate_galaxy{G <: AbstractFloat, T <: Integer}(p::SFWParameters{G}, samplesize::T; x0::G = 10.0^-8, rate::Bool = false)
    function sfwopt(x::Array{G, 1})
        return -(x[1] * x[1] * x[3] * profile_density(x[1], x[2], x[3], p))
    end

    minProb = Optim.optimize(sfwopt, [0.1, 1.0, 1.0], method = :nelder_mead)
    fmax = max(0.0,  -1.1 * minProb.f_minimum)
    
    vmax0 = escape_velocity(x0, p)
    sampledValues = Array{G}(samplesize, 3)
    if rate
        tested = 0
    end
    accepted = 0
    while accepted < samplesize
        x = (p.rlim  / p.rs) * rand()
        vr = vmax0 * rand()
        vt = vmax0 * rand()
        u = rand()
        sampleDensity = x * x * vt * profile_density(x, vr, vt, p)
        if rate
            tested += 1
        end
        if sampleDensity >= (u * fmax)
            accepted += 1
            sampledValues[accepted, :] = [x * p.rs , vr, vt]
        end
    end
    if rate
        return sampledValues, accepted / tested
    end
    return sampledVaues
end