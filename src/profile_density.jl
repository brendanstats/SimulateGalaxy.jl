"""
Compute the dark matter density profile at the point (x, vr, vt).  Where x is the radius divided by the scale radius (r / rs), vr is the radial velocity, and vt is the tangential velocity.
"""
function profile_density(x::G, vr::G, vt::G, p::NFWParameters{G}; N::G = 1000.0, αt::G = 1.0) where G <: AbstractFloat
    if x < 0.0
        return 0.0
    end

    #Compute Energy
    E = 0.5 * (vt * vt + vr * vr) + gravitational_potential(x, p)
    if E < p.Φlim && E >= 0.0
        hE = N * E^p.a * (E^p.q + p.adjEc^p.q)^(p.d/p.q) * (p.Φlim - E)^p.e
    else
        return 0.0
    end

    #Compute angular momentum
    J = abs(x * p.rs * vt) #vt = v * sin(theta)
    if p.b <= 0.0
        #gJ = 1.0 / (1.0 + (J / p.adjJb)^(-p.b))
        gJ = (1.0 + (J / p.adjJb)^(-p.b / αt))^-αt
    else
        #gJ = 1.0 + (J / p.adjJb)^p.b
        gJ = (1.0 + (J / p.adjJb)^(p.b / αt))^αt
    end

    return hE * gJ
end

function log_profile_density(x::G, vr::G, vt::G, p::NFWParameters{G}; N::G = 1000.0, αt::G = 1.0) where G <: AbstractFloat
    if x < 0.0
        return -Inf
    end

    #Compute Energy
    E = 0.5 * (vt * vt + vr * vr) + gravitational_potential(x, p)
    lE = log(E)
    if E < p.Φlim && E >= 0.0
        lhE = log(N) + p.a * lE + (p.d/p.q) * log(E^p.q + p.adjEc^p.q) + p.e * log(p.Φlim - E)
    else
        return -Inf
    end

    #Compute angular momentum
    J = abs(x * p.rs * vt) #vt = v * sin(theta)
    if p.b <= 0.0
        #gJ = 1.0 / (1.0 + (J / p.adjJb)^(-p.b))
        lgJ = -αt * log(1.0 + (J / p.adjJb)^(-p.b / αt))
    else
        #gJ = 1.0 + (J / p.adjJb)^p.b
        lgJ = αt * log(1.0 + (J / p.adjJb)^(p.b / αt))
    end

    return lhE + lgJ
end

function profile_density(x::G, vr::G, vt::G, p::SFWParameters{G}; N::G = 1000.0, αt::G = 1.0) where G <: AbstractFloat
    if x < 0.0
        return 0.0
    end

    #Compute Energy
    E = 0.5 * (vt * vt + vr * vr) + gravitational_potential(x, p)
    if E < p.Φlim && E >= 0.0
        hE = N * E^p.a * (E^p.q + p.adjEc^p.q)^(p.d/p.q) * (p.Φlim - E)^p.e
    else
        return 0.0
    end

    #Compute angular momentum
    J = abs(x * p.rs * vt) #vt = v * sin(theta)
    if p.b <= 0.0
        #gJ = 1.0 / (1.0 + (J / p.adjJb)^(-p.b))
        gJ = (1.0 + (J / p.adjJb)^(-p.b / αt))^-αt
    else
        #gJ = 1.0 + (J / p.adjJb)^p.b
        gJ = (1.0 + (J / p.adjJb)^(p.b / αt))^αt
    end

    return hE * gJ
end

function log_profile_density(x::G, vr::G, vt::G, p::SFWParameters{G}; N::G = 1000.0, αt::G = 1.0) where G <: AbstractFloat
    if x < 0.0
        return -Inf
    end

    #Compute Energy
    E = 0.5 * (vt * vt + vr * vr) + gravitational_potential(x, p)
    lE = log(E)
    if E < p.Φlim && E >= 0.0
        lhE = log(N) + p.a * lE + (p.d/p.q) * log(E^p.q + p.adjEc^p.q) + p.e * log(p.Φlim - E)
    else
        return -Inf
    end

    #Compute angular momentum
    J = abs(x * p.rs * vt) #vt = v * sin(theta)
    if p.b <= 0.0
        lgJ = -αt * log(1.0 + (J / p.adjJb)^(-p.b / αt))
    else
        lgJ = αt * log(1.0 + (J / p.adjJb)^(-p.b / αt))
    end

    return lhE + lgJ
end
