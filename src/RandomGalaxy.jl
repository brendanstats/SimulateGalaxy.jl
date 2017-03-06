abstract RandomGalaxy{G <: AbstractFloat, T <: Integer}

type SphericalGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    vr::Array{G, 1}
    vt::Array{G, 1}
    nobs::T
    #SphericalGalaxy(r, vr, vt, nobs) = (length(r) == nobs && length(vr) == nobs && length(vt) == nobs) ? error("length of r, vr, and vt must equal nobs") : new(r, vr, vt, nobs)
end

SphericalGalaxy{G <: AbstractFloat}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}) = SphericalGalaxy(r, vr, vt, length(r)) 

type EuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    z::Array{G, 1}
    vx::Array{G, 1}
    vy::Array{G, 1}
    vz::Array{G, 1}
    nobs::T
end

EuclideanGalaxy{G <: AbstractFloat}(x::Array{G, 1}, y::Array{G, 1}, z::Array{G, 1}, vx::Array{G, 1}, vy::Array{G, 1}, vz::Array{G, 1}) = EuclideanGalaxyGalaxy(x, y, z, vx, vy, vz, length(x)) 

type PartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    vz::Array{G, 1}
    nobs::T
end

PartialEuclideanGalaxy{G <: AbstractFloat}(x::Array{G, 1}, y::Array{G, 1}, vz::Array{G, 1}) = EuclideanGalaxyGalaxy(x, y, vz, length(x)) 


"""
Generate a random EuclideanGalaxy star from  a SphericalGalaxy star
"""
function sample_euclidean{G <: AbstractFloat}(r::G, vr::G, vt::G)
    theta = acos(1.0 - 2.0 * rand())
    phi = 2.0 * pi * rand()

    velocity = sqrt(vr^2 + vt^2)
    velSign = sign(rand() - 0.5)
    velPhi = 2.0 * pi * rand()

    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)

    vz2 = velSign * vr
    vx2 = vt * cos(velPhi)
    vy2 = vt * sin(velPhi)

    vx = cos(theta) * cos(phi) * vx2 - sin(phi) * vy2 + sin(theta) * cos(phi) * vz2
    vy = cos(theta) * sin(phi) * vx2 + cos(phi) * vy2 + sin(theta) * sin(phi) * vz2
    vz = -sin(theta) * vx2 + cos(theta) * vz2
    return x, y, z, vx, vy, vz
end

"""
Generate a EuclideanGalaxy from a SphericalGalaxy
"""
function sample_euclidean{G <: AbstractFloat}(sg::SphericalGalaxy{G})
    x = Array{G}(sg.nobs)
    y = Array{G}(sg.nobs)
    z = Array{G}(sg.nobs)
    vx = Array{G}(sg.nobs)
    vy = Array{G}(sg.nobs)
    vz = Array{G}(sg.nobs)

    for ii in 1:sg.nobs
        x[ii], y[ii], z[ii], vx[ii], vy[ii], vz[ii] = sample_euclidean(sg.r[ii], sg.vr[ii], sg.vt[ii])
    end
    return EuclideanGalaxy(x, y, z, vx, vy, vz, sg.nobs)
end

"""
Generate a random PartialEuclideanGalaxy star from  a SphericalGalaxy star
"""
function sample_partial_euclidean{G <: AbstractFloat}(r::G, vr::G, vt::G)
    theta = acos(1.0 - 2.0 * rand())
    phi = 2.0 * pi * rand()

    velocity = sqrt(vr^2 + vt^2)
    velSign = sign(rand() - 0.5)
    velPhi = 2.0 * pi * rand()

    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)

    vz2 = velSign * vr
    vx2 = vt * cos(velPhi)

    vz = -sin(theta) * vx2 + cos(theta) * vz2
    return x, y, vz
end

"""
Generate a PartialEuclideanGalaxy from a SphericalGalaxy
"""
function sample_partial_euclidean{G <: AbstractFloat}(sg::SphericalGalaxy{G})
    x = Array{G}(sg.nobs)
    y = Array{G}(sg.nobs)
    vz = Array{G}(sg.nobs)

    for ii in 1:sg.nobs
        x[ii], y[ii], vz[ii] = sample_euclidean(sg.r[ii], sg.vr[ii], sg.vt[ii])
    end
    return PartialEuclideanGalaxy(x, y, vz, sg.nobs)
end
