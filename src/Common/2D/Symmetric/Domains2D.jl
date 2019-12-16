function LatticeDomain(
            boundary::Boundary{1},
            lattice::Lattice{1},
            n::Number = 1;
            type::Symbol = :generic,
            name::Symbol = :anonymous,
            fit::Bool = true
            )

    bnd = boundary
    lat = lattice
    if fit
        nx = ceil(Int,bnd.shape.a/lat.dx)
        dx = bnd.shape.a/nx
        lat = Lattice(dx; origin=bnd.shape.start+dx/2)
    end

    i1 = latticeindex(lat,bnd.shape.start)
    i2 = latticeindex(lat,bnd.shape.stop)
    imin, imax = floor(Int,min(i1,i2)), floor(Int,max(i1,i2))
    x = lat[imin:imax]
    inds_keep = bnd.shape.(x)
    x = x[inds_keep]
    indices = map(CartesianIndex,(imin:imax)[inds_keep])
    interior = trues(size(x))
    surface = falses(size(x))
    surface[1] = surface[end] = true
    bulk = interior .& .!surface
    corner = falses(size(x))
    nnm = (vcat(0,1:length(x)-1),)
    nnp = (vcat(2:length(x),0),)

    return LatticeDomain(boundary, lat, n, type, name, x, indices, interior,
                            bulk, surface, corner, nnm, nnp)
end

LatticeDomain(lattice::Lattice{1}, boundary::Boundary{1}, args...; kwargs...) =
    Lattice(boundary, lattice, args...; kwargs...)
