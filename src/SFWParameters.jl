type SFWParameters{G <: AbstractFloat}
    a::G
    d::G
    e::G
    Ec::G
    rlim::G
    b::G
    q::G
    Jb::G
    ρs::G
    rs::G
    α::G
    β::G
    γ::G

    #secondary parameters
    Φs::G
    xlim::G
    Φ0::G
    Φlim::G

    adjEc::G
    adjJb::G
end

function SFWParameters{G <: AbstractFloat}(a::G, d::G, e::G, Ec::G, rlim::G, b::G, q::G, Jb::G, ρs::G, rs::G, α::G, β::G, γ::G; x0 = 10.0^-8)

    Φs = rhos * rs ^ 2
    xlim = rlim / rs
    adjEc = Φs * Ec
    adjJb = Jb * rs * sqrt(Φs)

    p = SFWParameters(a, d, e, Ec, rlim, b, q, Jb, ρs, rs, α, β, γ, Φs, xlim, 0.0, 0.0, adjEc, adjJb)

    p.Φ0 = gravitational_potential(x0, p)
    p.Φlim = gravitational_potential(xlim, p)
    return p
end
