# SimulateGalaxy

[![Build Status](https://travis-ci.org/brendanstats/SimulateGalaxy.jl.svg?branch=master)](https://travis-ci.org/brendanstats/SimulateGalaxy.jl)

[![Coverage Status](https://coveralls.io/repos/brendanstats/SimulateGalaxy.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/brendanstats/SimulateGalaxy.jl?branch=master)

[![codecov.io](http://codecov.io/github/brendanstats/SimulateGalaxy.jl/coverage.svg?branch=master)](http://codecov.io/github/brendanstats/SimulateGalaxy.jl?branch=master)

###Density Profiles

######Navarro–Frenk–White (NFW) profile
A type containing the parameters to define a NFW profile is included within the type `NFWParameters`

######Generalized Density profile
Support is also provided for a more general density profile as descrbed in

Geringer-Sameth, Alex, Savvas M. Koushiappas, and Matthew Walker. "Dwarf galaxy annihilation and decay emission profiles for dark matter experiments." The Astrophysical Journal 801.2 (2015): 74.
http://stacks.iop.org/0004-637X/801/i=2/a=74

The appropriate type for such a profile is `SFWParameters`

###Random Galaxy
The `RandomGalaxy{G <: AbstractFloat, T <: Integer}` is abstract type which governs output of
###Rejection Sampler