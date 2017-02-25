"""
Transform the results of GenGalaxy to Euclidean coordinates.
"""
function converttoeuclidean{G <: AbstractFloat}(coordinate::Array{G,1})
    radius = coordinate[1]
    theta = acos(1.0 - 2.0 * rand())
    phi = 2.0 * pi * rand()

    velocity = sqrt(coordinate[2]^2 + coordinate[3]^2)
    velSign = sign(rand() - 0.5)
    velPhi = 2.0 * pi * rand()

    x = radius * sin(theta) * cos(phi)
    y = radius * sin(theta) * sin(phi)
    z = radius * cos(theta)

    vz2 = velSign * coordinate[2]
    vx2 = coordinate[3] * cos(velPhi)
    vy2 = coordinate[3] * sin(velPhi)

    vx = cos(theta) * cos(phi) * vx2 - sin(phi) * vy2 + sin(theta) * cos(phi) * vz2
    vy = cos(theta) * sin(phi) * vx2 + cos(phi) * vy2 + sin(theta) * sin(phi) * vz2
    vz = -sin(theta) * vx2 + cos(theta) * vz2
    return x, y, z, vx, vy, vz
end
