module SimulateGalaxy

# package code goes here
import GSL.hypergeom
import QuadGK.quadgk
import Optim, StatsBase

export NFWParameters
export SFWParameters
export RandomGalaxy,
    SphericalGalaxy,
    MetallicSphericalGalaxy,
    PlaneGalaxy,
    MetallicPlaneGalaxy,
    EuclideanGalaxy,
    MetallicEuclideanGalaxy,
    PartialEuclideanGalaxy,
    MetallicPartialEuclideanGalaxy,
    sample_plane,
    sample_euclidean,
    sample_partial_euclidean
export hypergeom_2F1_expanded,
    density_profile,
    gravitational_potential,
    gravitational_potential_integral,
    gravitational_potential0,
    gravitational_potential0_robust
export escape_velocity
export profile_density
export simulate_galaxy

include("NFWParameters.jl")
include("SFWParameters.jl")
include("RandomGalaxy.jl")
include("gravitational_potential.jl")
include("escape_velocity.jl")
include("profile_density.jl")
include("simulate_galaxy.jl")

end # module
