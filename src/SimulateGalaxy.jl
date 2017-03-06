module SimulateGalaxy

# package code goes here
import GSL.hypergeom
import Optim

export NFWParameters
export SFWParameters
export RandomGalaxy,
    SphericalGalaxy,
    EuclideanGalaxy,
    PartialEuclideanGalaxy,
    sample_euclidean,
    sample_partial_euclidean
export hypergeom_2F1_expanded,
    density_profile,
    gravitational_potential,
    gravitational_potential_integral,
    gravitational_potential0
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
