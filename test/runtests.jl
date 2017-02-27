import SimulateGalaxy
using Base.Test
include("galaxy_generation.jl")

#Test NFW Definition
nfwMR = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5]
nfwMP = [2.4, -7.9, 1.1, 0.17, 3.0, 0.0, 8.2, 0.086, 21.0, 1.5]
nfwpMR = SimulateGalaxy.NFWParameters(nfwMR...)
nfwpMP = SimulateGalaxy.NFWParameters(nfwMP...)

@testset "Computed NFW Parameters" begin
    @testset "MR" begin
        @test nfwpMR.rs == nfwMR[10] / 2.16
        @test nfwpMR.Φs == (nfwMR[9] / 0.465) ^ 2
        @test nfwpMR.xlim == nfwMR[5] / nfwMR[10] * 2.16
        @test nfwpMR.Φlim == (nfwMR[9] / 0.465) ^ 2 * GravPotential(nfwMR[5] / nfwMR[10] * 2.16)
        @test nfwpMR.adjEc == nfwMR[4] * (nfwMR[9] / 0.465) ^ 2
        @test_approx_eq nfwpMR.adjJb nfwMR[8] * nfwMR[10] / 2.16 * sqrt((nfwMR[9] / 0.465) ^ 2)
    end
    @testset "MP" begin
        @test nfwpMP.rs == nfwMP[10] / 2.16
        @test nfwpMP.Φs == (nfwMP[9] / 0.465) ^ 2
        @test nfwpMP.xlim == nfwMP[5] / nfwMP[10] * 2.16
        @test nfwpMP.Φlim == (nfwMP[9] / 0.465) ^ 2 * GravPotential(nfwMP[5] / nfwMP[10] * 2.16)
        @test nfwpMP.adjEc == nfwMP[4] * (nfwMP[9] / 0.465) ^ 2
        @test_approx_eq nfwpMP.adjJb nfwMP[8] * nfwMP[10] / 2.16 * sqrt((nfwMP[9] / 0.465) ^ 2)
    end
end

#Test SFW Definition
sfwMR = [2.0, -5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 4229.2, .69444, 1.2, 3.05, 1.1]
sfwMP = [2.4, -7.9, 1.1, 0.17, 3.0, 0.0, 8.2, 0.0, 4229.2, .69444, 1.2, 3.05, 1.1]
sfwpMR = SimulateGalaxy.SFWParameters(sfwMR...)
sfwpMP = SimulateGalaxy.SFWParameters(sfwMP...)

@testset "Computed NFW Parameters" begin
    @testset "MR" begin
        @test sfwpMR.Φs == sfwMR[9] * sfwMR[10]^2
        @test sfwpMR.xlim == sfwMR[5] / sfwMR[10]
        @test sfwpMR.Φ0 == genphi0(10.0^-8, sfwMR[11], sfwMR[12], sfwMR[13], sfwMR[9], sfwMR[10])
        @test sfwpMR.Φlim == genphi(sfwMR[5] / sfwMR[10], sfwMR[11], sfwMR[12], sfwMR[13], sfwMR[9], sfwMR[10])
        @test sfwpMR.adjEc == sfwMR[4] * sfwMR[9] * sfwMR[10]^2
        @test_approx_eq sfwpMR.adjJb sfwMR[8] * sqrt(sfwMR[9]) * sfwMR[10]^2
    end
    @testset "MP" begin
        @test sfwpMP.Φs == sfwMP[9] * sfwMP[10]^2
        @test sfwpMP.xlim == sfwMP[5] / sfwMP[10]
        @test sfwpMP.Φ0 == genphi0(10.0^-8, sfwMP[11], sfwMP[12], sfwMP[13], sfwMP[9], sfwMP[10])
        @test sfwpMP.Φlim == genphi(sfwMP[5] / sfwMP[10], sfwMP[11], sfwMP[12], sfwMP[13], sfwMP[9], sfwMP[10])
        @test sfwpMP.adjEc == sfwMP[4] * sfwMP[9] * sfwMP[10]^2
        @test_approx_eq sfwpMP.adjJb sfwMP[8] * sqrt(sfwMP[9]) * sfwMP[10]^2
    end
end
    

#Test gravitational_potential
n = 5
xMR = rand(n) * nfwpMR.rlim / nfwpMR.rs
xMP = rand(n) * nfwpMP.rlim / nfwpMP.rs
@testset "gravitational_potential Tests" begin
    @testset "MR" begin
        for x in xMR
            @test SimulateGalaxy.gravitational_potential(x) == GravPotential(x)
            @test SimulateGalaxy.gravitational_potential(x, nfwpMR) == GravPotential(x) * (nfwMR[9] / 0.465) ^ 2
            @test SimulateGalaxy.gravitational_potential(x, sfwpMR) == genphi(x, sfwMR[11], sfwMR[12], sfwMR[13], sfwMR[9], sfwMR[10])
        end
    end
    @testset "MP" begin
        for x in xMP
            @test SimulateGalaxy.gravitational_potential(x) == GravPotential(x)
            @test SimulateGalaxy.gravitational_potential(x, nfwpMP) == GravPotential(x) * (nfwMP[9] / 0.465) ^ 2
            @test SimulateGalaxy.gravitational_potential(x, sfwpMP) == genphi(x, sfwMP[11], sfwMP[12], sfwMP[13], sfwMP[9], sfwMP[10])
        end
    end
end

#Test escape_velocity
@testset "escape_velocity Tests" begin
    @testset "MR" begin
        for x in xMR
            @test_approx_eq SimulateGalaxy.escape_velocity(x, nfwpMR) EscapeVelocity(x, nfwMR[5], nfwMR[10] / 2.16, (nfwMR[9] / 0.465) ^ 2)
            @test SimulateGalaxy.escape_velocity(x, sfwpMR) == escape_velocity(x, sfwMR[5], sfwMR[9], sfwMR[10], sfwMR[11], sfwMR[12], sfwMR[13])
        end
    end
    @testset "MP" begin
        for x in xMP
            @test_approx_eq SimulateGalaxy.escape_velocity(x, nfwpMP) EscapeVelocity(x, nfwMP[5], nfwMP[10] / 2.16, (nfwMP[9] / 0.465) ^ 2)
            @test SimulateGalaxy.escape_velocity(x, sfwpMP) == escape_velocity(x, sfwMP[5], sfwMP[9], sfwMP[10], sfwMP[11], sfwMP[12], sfwMP[13])
        end
    end
end

#Test NFW profile_density
vrMR = rand(n) * nfwpMR.vmax
vtMR = rand(n) * nfwpMR.vmax
vrMP = rand(n) * nfwpMP.vmax
vtMP = rand(n) * nfwpMP.vmax

@testset "NFW profile_density tests" begin
    @testset "MR" begin
        for (x, vr, vt) in zip(xMR, vrMR, vtMR)
            @test_approx_eq SimulateGalaxy.profile_density(x, vr, vt, nfwpMR) NFWDensity(x, vr, vt, nfwMR)
        end
    end
    @testset "MP" begin
        for (x, vr, vt) in zip(xMP, vrMP, vtMP)
            @test SimulateGalaxy.profile_density(x, vr, vt, nfwpMP) == NFWDensity(x, vr, vt, nfwMP)
        end
    end
end

#Test SFW profile_density
@testset "SFW profile_density tests" begin
    @testset "MR" begin
        for (x, vr, vt) in zip(xMR, vrMR, vtMR)
            @test_approx_eq SimulateGalaxy.profile_density(x, vr, vt, sfwpMR) SFW_density(x, vr, vt, sfwMR)
        end
    end
    @testset "MP" begin
        for (x, vr, vt) in zip(xMP, vrMP, vtMP)
            @test SimulateGalaxy.profile_density(x, vr, vt, sfwpMP) == SFW_density(x, vr, vt, sfwMP)
        end
    end
end

#Test simulate_galaxy NFW
@time SimulateGalaxy.simulate_galaxy(nfwpMR, 400, rate = true)
@time SimulateGalaxy.simulate_galaxy(nfwpMP, 400, rate = true)
@time SimulateGalaxy.simulate_galaxy(sfwpMR, 400, rate = true)
@time SimulateGalaxy.simulate_galaxy(sfwpMP, 400, rate = true)
