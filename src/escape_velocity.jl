"""
Calculate escape velocity for a radius scaled by the
characteristic radius, i.e. x = r / rs
"""
function escape_velocity(x::G, p::NFWParameters) where G <: AbstractFloat
    return sqrt(2.0  * (p.Φlim - gravitational_potential(x, p)))
end

function escape_velocity(x::G, p::SFWParameters) where G <: AbstractFloat
    return sqrt(2.0  * (p.Φlim - gravitational_potential(x, p)))
end
