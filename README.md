# SimulateGalaxy

[![Build Status](https://travis-ci.org/brendanstats/SimulateGalaxy.jl.svg?branch=master)](https://travis-ci.org/brendanstats/SimulateGalaxy.jl)

[![Coverage Status](https://coveralls.io/repos/brendanstats/SimulateGalaxy.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/brendanstats/SimulateGalaxy.jl?branch=master)

[![codecov.io](http://codecov.io/github/brendanstats/SimulateGalaxy.jl/coverage.svg?branch=master)](http://codecov.io/github/brendanstats/SimulateGalaxy.jl?branch=master)

### Density Profiles

###### Navarro–Frenk–White (NFW) profile
A type containing the parameters to define a NFW profile is included within the type `NFWParameters`

###### Generalized Density profile
Support is also provided for a more general density profile as descrbed in

Geringer-Sameth, Alex, Savvas M. Koushiappas, and Matthew Walker. "Dwarf galaxy annihilation and decay emission profiles for dark matter experiments." The Astrophysical Journal 801.2 (2015): 74.
http://stacks.iop.org/0004-637X/801/i=2/a=74

The appropriate type for such a profile is `SFWParameters`

### RandomGalaxy
The `RandomGalaxy{G <: AbstractFloat, T <: Integer}` is abstract type which encompasses the different types of output that can be generated.  The subtypes which are generated and the fields they contain are as follows:
* `SphericalGalaxy{G <: AbstractFloat, T <: Integer}`
  + r
  + vr
  + vt
  + nobs
* `MetallicSphericalGalaxy{G <: AbstractFloat, T <: Integer}`
  + r
  + vr
  + vt
  + m
  + nobs
* `EuclideanGalaxy{G <: AbstractFloat, T <: Integer}`
  + x
  + y
  + z
  + vx
  + vy
  + vz
  + nobs
* `MetallicEuclideanGalaxy{G <: AbstractFloat, T <: Integer}`
  + x
  + y
  + z
  + vx
  + vy
  + vz
  + m
  + nobs
* `PartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer}`
  + x
  + y
  + vz
  + nobs
* `MetallicPartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer}`
  + x
  + y
  + vz
  + m
  + nobs

`SphericalGalaxy` and `MetallicSphericalGalaxy` denote a spherical coordinate space where only radius, radial velocity, and tangential velocity are specified.  `EuclideanGalaxy` and `MetallicEuclideanGalaxy` are embeded in Euclidean space for both position and veloicty.  Finally `PartialEuclideanGalaxy` and `MetallicPartialEuclideanGalaxy` denote imperfectly observed galaxies in Euclidean space containing only two coordinates for position and one of velocity.  In all cases the Metallic prefix indicates that a `m` field denoting metallicity is included.

### Rejection Sampler
A rejection sampler is used to draw from the density profile returning a radius, radial velocity, and tangential velocity.  If parameters defining the mean and standard deviation of a metallicity are included then a normally distributed metallicity observation is generated as well.  Simulation is done via the `simulate_galaxy` function.

### Examples
```{julia}
using SimulateGalaxy

#Define Density Profiles
nfw = NFWParameters(2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5)
sfw = SFWParameters(2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 4229.2, .69444, 1.2, 3.05, 1.1)

#Simulate Spherical Galaxies with 400 Stars
nfw_sph = simulate_galaxy(nfw, 400, rate = true)
sfw_sph = simulate_galaxy(sfw, 400, rate = true)

#Draw a random conversion to Euclidean coordinates based on the spherical coordinates
nfw_eu = sample_euclidean(nfw_sph)
sfw_eu = sample_euclidean(sfw_sph)

#Simulate Spherical Galaxies with 400 Stars and a metallicity component
nfw_sph = simulate_galaxy(nfw, -2.0, .003, 400, rate = true)
sfw_sph = simulate_galaxy(sfw, -2.0, .003, 400, rate = true)
```
