"""
Parameters for NFW Profile, includes useful quantities defined by parameters
"""
type NFWParameters{G <: AbstractFloat}
    a::G
    d::G
    e::G
    Ec::G
    rlim::G
    b::G
    q::G
    Jb::G
    vmax::G
    rmax::G

    #secondary parameters
    rs::G
    Φs::G
    xlim::G
    Φlim::G

    #adjusted parameter values
    adjEc::G
    adjJb::G
end

function NFWParameters(a::G, d::G, e::G, Ec::G, rlim::G, b::G, q::G, Jb::G, vmax::G, rmax::G) where G <: AbstractFloat

    rs = rmax / 2.16
    Φs = (vmax / 0.465)^2
    xlim = rlim / rs
    adjEc = Φs * Ec
    adjJb = Jb * rs * vmax / 0.465

    p = NFWParameters(a, d, e, Ec, rlim, b, q, Jb, vmax, rmax, rs, Φs, xlim, 0.0, adjEc, adjJb)
    p.Φlim = gravitational_potential(xlim, p)
    
    return p
end
