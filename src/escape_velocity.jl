"""
Calculate escape velocity for a radius scaled by the
characteristic radius, i.e. x = r / rs
"""
function escape_velocity{G <: AbstractFloat}(x::G, p::NFWParameters)
    return sqrt(2.0  * (p.Φlim - gravitational_potential(x, p)))
end

function escape_velocity{G <: AbstractFloat}(x::G, p::SFWParameters)
    return sqrt(2.0  * (p.Φlim - gravitational_potential(x, p)))
end
