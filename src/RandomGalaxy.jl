abstract type RandomGalaxy{G <: AbstractFloat, T <: Integer} end

type SphericalGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    vr::Array{G, 1}
    vt::Array{G, 1}
    nobs::T
    #SphericalGalaxy(r, vr, vt, nobs) = (length(r) == nobs && length(vr) == nobs && length(vt) == nobs) ? error("length of r, vr, and vt must equal nobs") : new(r, vr, vt, nobs)
end

SphericalGalaxy(sg1::SphericalGalaxy{G, T}, sg2::SphericalGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = SphericalGalaxy([sg1.r; sg2.r], [sg1.vr; sg2.vr], [sg1.vt; sg2.vt], sg1.nobs + sg2.nobs)
SphericalGalaxy(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}) where G <: AbstractFloat = SphericalGalaxy(r, vr, vt, length(r))

type MetallicSphericalGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    vr::Array{G, 1}
    vt::Array{G, 1}
    m::Array{G, 1}
    nobs::T
    #SphericalGalaxy(r, vr, vt, nobs) = (length(r) == nobs && length(vr) == nobs && length(vt) == nobs) ? error("length of r, vr, and vt must equal nobs") : new(r, vr, vt, nobs)
end

MetallicSphericalGalaxy(msg1::MetallicSphericalGalaxy{G, T}, msg2::MetallicSphericalGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = MetallicSphericalGalaxy([msg1.r; msg2.r], [msg1.vr; msg2.vr], [msg1.vt; msg2.vt], [msg1.m; msg2.m], msg1.nobs + msg2.nobs)
MetallicSphericalGalaxy(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}, m::Array{G, 1}) where G <: AbstractFloat = MetallicSphericalGalaxy(r, vr, vt, m, length(r))
MetallicSphericalGalaxy(sg::SphericalGalaxy{G, T}, m::Array{G,1}) where {G <: AbstractFloat, T <: Integer} = MetallicSphericalGalaxy(sg.r, sg.vr, sg.vt, m, sg.nobs) 

type PlaneGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    v::Array{G, 1}
    nobs::T
end

PlaneGalaxy(dg1::PlaneGalaxy{G, T}, dg2::PlaneGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = PlaneGalaxy([dg1.r; dg2.r], [dg1.v; dg2.v], dg1.nobs + dg2.nobs)
PlaneGalaxy(r::Array{G, 1}, v::Array{G, 1}) where G <: AbstractFloat = PlaneGalaxy(r, v, length(r))


type MetallicPlaneGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    r::Array{G, 1}
    v::Array{G, 1}
    m::Array{G, 1}
    nobs::T
end

MetallicPlaneGalaxy(mdg1::PlaneGalaxy{G, T}, mdg2::PlaneGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = MetallicPlaneGalaxy([mdg1.r; mdg2.r], [mdg1.v; mdg2.v], [mdg1.m; mdg2.m], dg1.nobs + dg2.nobs)
MetallicPlaneGalaxy(r::Array{G, 1}, v::Array{G, 1}, m::Array{G, 1}) where G <: AbstractFloat = PlaneGalaxy(r, v, m, length(r))
MetallicPlaneGalaxy(mdg::PlaneGalaxy{G, T}, m::Array{G,1}) where {G <: AbstractFloat, T <: Integer} = MetallicPlaneGalaxy(mdg.r, mdg.v, m, mdg.nobs) 

type EuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    z::Array{G, 1}
    vx::Array{G, 1}
    vy::Array{G, 1}
    vz::Array{G, 1}
    nobs::T
end

EuclideanGalaxy(eg1::EuclideanGalaxy{G, T}, eg2::EuclideanGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = Euclideanalaxy([eg1.x; eg2.x], [eg1.y; eg2.y], [eg1.z; eg2.z], [eg1.vx; eg2.vx], [eg1.vy; eg2.vy], [eg1.vz; eg2.vz], eg1.nobs + eg2.nobs)
EuclideanGalaxy(x::Array{G, 1}, y::Array{G, 1}, z::Array{G, 1}, vx::Array{G, 1}, vy::Array{G, 1}, vz::Array{G, 1}) where G <: AbstractFloat = EuclideanGalaxy(x, y, z, vx, vy, vz, length(x)) 

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

MetallicEuclideanGalaxy(meg1::MetallicEuclideanGalaxy{G, T}, meg2::MetallicEuclideanGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = MetallicEuclideanGalaxy([meg1.x; meg2.x], [meg1.y; meg2.y],[meg1.z; meg2.z], [meg1.vx; meg2.vx], [meg1.vy; meg2.vy], [meg1.vz; meg2.vz], [meg1.m; meg2.m], meg1.nobs + meg2.nobs)
MetallicEuclideanGalaxy(x::Array{G, 1}, y::Array{G, 1}, z::Array{G, 1}, vx::Array{G, 1}, vy::Array{G, 1}, vz::Array{G, 1}, m::Array{G, 1}) where G <: AbstractFloat= MetallicEuclideanGalaxy(x, y, z, vx, vy, vz, length(x))
MetallicEuclideanGalaxy{G <: AbstractFloat, T <: Integer}(eg::EuclideanGalaxy{G, T}, m::Array{G, 1}) = MetallicEuclideanGalaxy(eg.x, eg.y, eg.z, eg.vx, eg.vy, eg.vz, m, eg.nobs) 

type PartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    vz::Array{G, 1}
    nobs::T
end

PartialEuclideanGalaxy(peg1::PartialEuclideanGalaxy{G, T}, peg2::PartialEuclideanGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = PartialEuclideanGalaxy([peg1.x; peg2.x], [peg1.y; peg2.y], [peg1.vz; peg2.vz], peg1.nobs + peg2.nobs)
PartialEuclideanGalaxy(x::Array{G, 1}, y::Array{G, 1}, vz::Array{G, 1}) where G <: AbstractFloat = PartialEuclideanGalaxy(x, y, vz, length(x))
PartialEuclideanGalaxy(eg::EuclideanGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = PartialEuclideanGalaxy(eg.x, eg.y, eg.vz, eg.nobs)

type MetallicPartialEuclideanGalaxy{G <: AbstractFloat, T <: Integer} <: RandomGalaxy{G, T}
    x::Array{G, 1}
    y::Array{G, 1}
    vz::Array{G, 1}
    m::Array{G, 1}
    nobs::T
end

MetallicPartialEuclideanGalaxy(mpeg1::MetallicPartialEuclideanGalaxy{G, T}, mpeg2::MetallicPartialEuclideanGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer} = MetallicPartialEuclideanGalaxy([mpeg1.x; mpeg2.x], [mpeg1.y; mpeg2.y], [mpeg1.vz; mpeg2.vz],  [mpeg1.m; mpeg2.m], mpeg1.nobs + mpeg2.nobs)
MetallicPartialEuclideanGalaxy(x::Array{G, 1}, y::Array{G, 1}, vz::Array{G, 1}, m::Array{G, 1}) where G <: AbstractFloat = MetallicPartialEuclideanGalaxy(x, y, vz, m, length(x))
MetallicPartialEuclideanGalaxy(peg::PartialEuclideanGalaxy{G, T}, m::Array{G, 1}) where {G <: AbstractFloat, T <: Integer} = MetallicPartialEuclideanGalaxy(peg.x, peg.y, peg.vz, m, peg.nobs)
MetallicPartialEuclideanGalaxy(eg::EuclideanGalaxy{G, T}, m::Array{G, 1}) where {G <: AbstractFloat, T <: Integer} = PartialEuclideanGalaxy(eg.x, eg.y, eg.vz, m, eg.nobs)

"""
Generate a random PlaneGalaxy star from  a SphericalGalaxy star
"""
function sample_plane(r::G, vr::G, vt::G) where G <: AbstractFloat
    u = rand()
    s = 2.0 * sqrt(u - u * u)

    velPhi = 2pi * rand()
    velSign = sign(rand() - 0.5)
    
    return r * s, -s * cos(velPhi) * vt + (1.0 - 2.0 * u) * velSign * vr
end

function sample_plane(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}, nobs::T) where {G <: AbstractFloat, T <: Integer}
    rp = Array{G}(nobs)
    vp = Array{G}(nobs)
    for (ii, (ri, vri, vti)) in enumerate(zip(r, vr, vt))
        rp[ii], vp[ii] = sample_plane(ri, vri, vti)
    end
    return rp, vp
end

function sample_plane(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}) where G <: AbstractFloat
    return sample_plane(r, vr, vt, length(r))
end

function sample_plane(sg::SphericalGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer}
    return PlaneGalaxy(sample_plane(sg.r, sg.vr, sg.vt, sg.nobs)..., sg.nobs)
end

function sample_plane(msg::MetallicSphericalGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer}
    return MetallicPlaneGalaxy(sample_plane(msg.r, msg.vr, msg.vt, msg.nobs)..., msg.m, msg.nobs)
end


"""
Generate a random EuclideanGalaxy star from  a SphericalGalaxy star
"""
function sample_euclidean(r::G, vr::G, vt::G) where G <: AbstractFloat
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

function sample_euclidean(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}, nobs::T) where {G <: AbstractFloat, T <: Integer}
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

function sample_euclidean(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}) where G <: AbstractFloat
    return sample_euclidean(r, vr, vt, length(r))
end

"""
Generate a EuclideanGalaxy from a SphericalGalaxy
"""
function sample_euclidean(sg::SphericalGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer}
    return EuclideanGalaxy(sample_euclidean(sg.r, sg.vr, sg.vt, sg.nobs)..., sg.nobs)
end

function sample_euclidean(msg::MetallicSphericalGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer}
    return MetallicEuclideanGalaxy(sample_euclidean(msg.r, msg.vr, msg.vt, msg.nobs)..., msg.m, msg.nobs)
end

"""
Generate a random PartialEuclideanGalaxy star from  a SphericalGalaxy star
"""
function sample_partial_euclidean(r::G, vr::G, vt::G) where G <: AbstractFloat
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

function sample_partial_euclidean(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}, nobs::T) where {G <: AbstractFloat, T <: Integer}
    x = Array{G}(nobs)
    y = Array{G}(nobs)
    vz = Array{G}(nobs)
    for (ii, (ri, vri, vti)) in enumerate(zip(r, vr, vt))
        x[ii], y[ii], vz[ii] = sample_partial_euclidean(ri, vri, vti)
    end
    return x, y, vz
end

function sample_partial_euclidean(r::Array{G, 1}, vr::Array{G, 1}, vt::Array{G, 1}) where G <: AbstractFloat
    return sample_euclidean(r, vr, vt, length(r))
end

"""
Generate a PartialEuclideanGalaxy from a SphericalGalaxy
"""
function sample_partial_euclidean(sg::SphericalGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer}
    return PartialEuclideanGalaxy(sample_partial_euclidean(sg.r, sg.vr, sg.vt, sg.nobs)..., sg.nobs)
end

function sample_partial_euclidean(msg::MetallicSphericalGalaxy{G, T}) where {G <: AbstractFloat, T <: Integer}
    return MetallicPartialEuclideanGalaxy(sample_partial_euclidean(msg.r, msg.vr, msg.vt, msg.nobs)..., msg.m, msg.nobs)
end
