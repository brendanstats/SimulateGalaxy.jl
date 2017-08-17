abstract type RandomGalaxy{G <: AbstractFloat, T <: Integer} end

type SphericalGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    vr::Array{G, 1}
    vt::Array{G, 1}
    nobs::T
    #SphericalGalaxy(r, vr, vt, nobs) = (length(r) == nobs && length(vr) == nobs && length(vt) == nobs) ? error("length of r, vr, and vt must equal nobs") : new(r, vr, vt, nobs)
end

SphericalGalaxy{G <: AbstractFloat, T <: Integer}(sg1::SphericalGalaxy{G, T}, sg2::SphericalGalaxy{G, T}) = SphericalGalaxy([sg1.r; sg2.r], [sg1.vr; sg2.vr], [sg1.vt; sg2.vt], sg1.nobs + sg2.nobs)
SphericalGalaxy{G <: AbstractFloat}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}) = SphericalGalaxy(r, vr, vt, length(r))

type MetallicSphericalGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    vr::Array{G, 1}
    vt::Array{G, 1}
    m::Array{G, 1}
    nobs::T
    #SphericalGalaxy(r, vr, vt, nobs) = (length(r) == nobs && length(vr) == nobs && length(vt) == nobs) ? error("length of r, vr, and vt must equal nobs") : new(r, vr, vt, nobs)
end

MetallicSphericalGalaxy{G <: AbstractFloat, T <: Integer}(msg1::MetallicSphericalGalaxy{G, T}, msg2::MetallicSphericalGalaxy{G, T}) = MetallicSphericalGalaxy([msg1.r; msg2.r], [msg1.vr; msg2.vr], [msg1.vt; msg2.vt], [msg1.m; msg2.m], msg1.nobs + msg2.nobs)
MetallicSphericalGalaxy{G <: AbstractFloat}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}, m::Array{G, 1}) = MetallicSphericalGalaxy(r, vr, vt, m, length(r))
MetallicSphericalGalaxy{G <: AbstractFloat, T <: Integer}(sg::SphericalGalaxy{G, T}, m::Array{G,1}) = MetallicSphericalGalaxy(sg.r, sg.vr, sg.vt, m, sg.nobs) 

type PlaneGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    v::Array{G, 1}
    nobs::T
end

PlaneGalaxy{G <: AbstractFloat, T <: Integer}(dg1::PlaneGalaxy{G, T}, dg2::PlaneGalaxy{G, T}) = PlaneGalaxy([dg1.r; dg2.r], [dg1.v; dg2.v], dg1.nobs + dg2.nobs)
PlaneGalaxy{G <: AbstractFloat}(r::Array{G, 1}, v::Array{G, 1}) = PlaneGalaxy(r, v, length(r))


type MetallicPlaneGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    v::Array{G, 1}
    m::Array{G, 1}
    nobs::T
end

MetallicPlaneGalaxy{G <: AbstractFloat, T <: Integer}(mdg1::PlaneGalaxy{G, T}, mdg2::PlaneGalaxy{G, T}) = MetallicPlaneGalaxy([mdg1.r; mdg2.r], [mdg1.v; mdg2.v], [mdg1.m; mdg2.m], dg1.nobs + dg2.nobs)
MetallicPlaneGalaxy{G <: AbstractFloat}(r::Array{G, 1}, v::Array{G, 1}, m::Array{G, 1}) = PlaneGalaxy(r, v, m, length(r))
MetallicPlaneGalaxy{G <: AbstractFloat, T <: Integer}(mdg::PlaneGalaxy{G, T}, m::Array{G,1}) = MetallicPlaneGalaxy(mdg.r, mdg.v, m, mdg.nobs) 

type EuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    z::Array{G, 1}
    vx::Array{G, 1}
    vy::Array{G, 1}
    vz::Array{G, 1}
    nobs::T
end

EuclideanGalaxy{G <: AbstractFloat, T <: Integer}(eg1::EuclideanGalaxy{G, T}, eg2::EuclideanGalaxy{G, T}) = Euclideanalaxy([eg1.x; eg2.x], [eg1.y; eg2.y], [eg1.z; eg2.z], [eg1.vx; eg2.vx], [eg1.vy; eg2.vy], [eg1.vz; eg2.vz], eg1.nobs + eg2.nobs)
EuclideanGalaxy{G <: AbstractFloat}(x::Array{G, 1}, y::Array{G, 1}, z::Array{G, 1}, vx::Array{G, 1}, vy::Array{G, 1}, vz::Array{G, 1}) = EuclideanGalaxy(x, y, z, vx, vy, vz, length(x)) 

type MetallicEuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    z::Array{G, 1}
    vx::Array{G, 1}
    vy::Array{G, 1}
    vz::Array{G, 1}
    m::Array{G, 1}
    nobs::T
end

MetallicEuclideanGalaxy{G <: AbstractFloat, T <: Integer}(meg1::MetallicEuclideanGalaxy{G, T}, meg2::MetallicEuclideanGalaxy{G, T}) = MetallicEuclideanGalaxy([meg1.x; meg2.x], [meg1.y; meg2.y],[meg1.z; meg2.z], [meg1.vx; meg2.vx], [meg1.vy; meg2.vy], [meg1.vz; meg2.vz], [meg1.m; meg2.m], meg1.nobs + meg2.nobs)
MetallicEuclideanGalaxy{G <: AbstractFloat}(x::Array{G, 1}, y::Array{G, 1}, z::Array{G, 1}, vx::Array{G, 1}, vy::Array{G, 1}, vz::Array{G, 1}, m::Array{G, 1}) = MetallicEuclideanGalaxy(x, y, z, vx, vy, vz, length(x))
MetallicEuclideanGalaxy{G <: AbstractFloat, T <: Integer}(eg::EuclideanGalaxy{G, T}, m::Array{G, 1}) = MetallicEuclideanGalaxy(eg.x, eg.y, eg.z, eg.vx, eg.vy, eg.vz, m, eg.nobs) 

type PartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    vz::Array{G, 1}
    nobs::T
end

PartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer}(peg1::PartialEuclideanGalaxy{G, T}, peg2::PartialEuclideanGalaxy{G, T}) = PartialEuclideanGalaxy([peg1.x; peg2.x], [peg1.y; peg2.y], [peg1.vz; peg2.vz], peg1.nobs + peg2.nobs)
PartialEuclideanGalaxy{G <: AbstractFloat}(x::Array{G, 1}, y::Array{G, 1}, vz::Array{G, 1}) = PartialEuclideanGalaxy(x, y, vz, length(x))
PartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer}(eg::EuclideanGalaxy{G, T}) = PartialEuclideanGalaxy(eg.x, eg.y, eg.vz, eg.nobs)

type MetallicPartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    vz::Array{G, 1}
    m::Array{G, 1}
    nobs::T
end

MetallicPartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer}(mpeg1::MetallicPartialEuclideanGalaxy{G, T}, mpeg2::MetallicPartialEuclideanGalaxy{G, T}) = MetallicPartialEuclideanGalaxy([mpeg1.x; mpeg2.x], [mpeg1.y; mpeg2.y], [mpeg1.vz; mpeg2.vz],  [mpeg1.m; mpeg2.m], mpeg1.nobs + mpeg2.nobs)
MetallicPartialEuclideanGalaxy{G <: AbstractFloat}(x::Array{G, 1}, y::Array{G, 1}, vz::Array{G, 1}, m::Array{G, 1}) = MetallicPartialEuclideanGalaxy(x, y, vz, m, length(x))
MetallicPartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer}(peg::PartialEuclideanGalaxy{G, T}, m::Array{G, 1}) = MetallicPartialEuclideanGalaxy(peg.x, peg.y, peg.vz, m, peg.nobs)
MetallicPartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer}(eg::EuclideanGalaxy{G, T}, m::Array{G, 1}) = PartialEuclideanGalaxy(eg.x, eg.y, eg.vz, m, eg.nobs)

"""
Generate a random PlaneGalaxy star from  a SphericalGalaxy star
"""
function sample_plane{G <: AbstractFloat}(r::G, vr::G, vt::G)
    u = rand()
    s = 2.0 * sqrt(u - u * u)

    velPhi = 2pi * rand()
    velSign = sign(rand() - 0.5)
    
    return r * s, -s * cos(velPhi) * vt + (1.0 - 2.0 * u) * velSign * vr
end

function sample_plane{G <: AbstractFloat, T <: Integer}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}, nobs::T)
    rp = Array{G}(nobs)
    vp = Array{G}(nobs)
    for (ii, (ri, vri, vti)) in enumerate(zip(r, vr, vt))
        rp[ii], vp[ii] = sample_plane(ri, vri, vti)
    end
    return rp, vp
end

function sample_plane{G <: AbstractFloat}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1})
    return sample_plane(r, vr, vt, length(r))
end

function sample_plane{G <: AbstractFloat, T <: Integer}(sg::SphericalGalaxy{G, T})
    return PlaneGalaxy(sample_plane(sg.r, sg.vr, sg.vt, sg.nobs)..., sg.nobs)
end

function sample_plane{G <: AbstractFloat, T <: Integer}(msg::MetallicSphericalGalaxy{G, T})
    return MetallicPlaneGalaxy(sample_plane(msg.r, msg.vr, msg.vt, msg.nobs)..., msg.m, msg.nobs)
end


"""
Generate a random EuclideanGalaxy star from  a SphericalGalaxy star
"""
function sample_euclidean{G <: AbstractFloat}(r::G, vr::G, vt::G)
    theta = acos(1.0 - 2.0 * rand())
    phi = 2.0 * pi * rand()

    #velocity = sqrt(vr^2 + vt^2)
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

function sample_euclidean{G <: AbstractFloat, T <: Integer}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}, nobs::T)
    x = Array{G}(nobs)
    y = Array{G}(nobs)
    z = Array{G}(nobs)
    vx = Array{G}(nobs)
    vy = Array{G}(nobs)
    vz = Array{G}(nobs)
    for (ii, (ri, vri, vti)) in enumerate(zip(r, vr, vt))
        x[ii], y[ii], z[ii], vx[ii], vy[ii], vz[ii] = sample_euclidean(ri, vri, vti)
    end
    return x, y, z, vx, vy, vz
end

function sample_euclidean{G <: AbstractFloat}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1})
    return sample_euclidean(r, vr, vt, length(r))
end

"""
Generate a EuclideanGalaxy from a SphericalGalaxy
"""
function sample_euclidean{G <: AbstractFloat, T <: Integer}(sg::SphericalGalaxy{G, T})
    return EuclideanGalaxy(sample_euclidean(sg.r, sg.vr, sg.vt, sg.nobs)..., sg.nobs)
end

function sample_euclidean{G <: AbstractFloat, T <: Integer}(msg::MetallicSphericalGalaxy{G, T})
    return MetallicEuclideanGalaxy(sample_euclidean(msg.r, msg.vr, msg.vt, msg.nobs)..., msg.m, msg.nobs)
end

"""
Generate a random PartialEuclideanGalaxy star from  a SphericalGalaxy star
"""
function sample_partial_euclidean{G <: AbstractFloat}(r::G, vr::G, vt::G)
    theta = acos(1.0 - 2.0 * rand())
    phi = 2pi * rand()

    velocity = sqrt(vr^2 + vt^2)
    velSign = sign(rand() - 0.5)
    velPhi = 2pi * rand()

    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)

    vz2 = velSign * vr
    vx2 = vt * cos(velPhi)

    vz = -sin(theta) * vx2 + cos(theta) * vz2
    return x, y, vz
end

function sample_partial_euclidean{G <: AbstractFloat, T <: Integer}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}, nobs::T)
    x = Array{G}(nobs)
    y = Array{G}(nobs)
    vz = Array{G}(nobs)
    for (ii, (ri, vri, vti)) in enumerate(zip(r, vr, vt))
        x[ii], y[ii], vz[ii] = sample_partial_euclidean(ri, vri, vti)
    end
    return x, y, vz
end

function sample_partial_euclidean{G <: AbstractFloat}(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1})
    return sample_euclidean(r, vr, vt, length(r))
end

"""
Generate a PartialEuclideanGalaxy from a SphericalGalaxy
"""
function sample_partial_euclidean{G <: AbstractFloat, T <: Integer}(sg::SphericalGalaxy{G, T})
    return PartialEuclideanGalaxy(sample_partial_euclidean(sg.r, sg.vr, sg.vt, sg.nobs)..., sg.nobs)
end

function sample_partial_euclidean{G <: AbstractFloat, T <: Integer}(msg::MetallicSphericalGalaxy{G, T})
    return MetallicPartialEuclideanGalaxy(sample_partial_euclidean(msg.r, msg.vr, msg.vt, msg.nobs)..., msg.m, msg.nobs)
end
