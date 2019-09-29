function Domain(
    boundary::Boundary{Interval},
    lattice::Lattice{1},
    dielectric::TDF = DielectricFunction(1),
    pump::TPF = PumpFunction(0);
    type::Symbol = :generic,
    name::Symbol = :anonymous,
    fit::Bool = false
    ) where
        {TDF<:DielectricFunction,
        TPF<:PumpFunction}

    bnd = boundary
    lat = lattice
    if fit
        nx = ceil(Int,bnd.shape.a/lat.dx)
        dx = bnd.shape.a/nx
        lat = Lattice(dx;origin=dx/2)
    end

    # if typeof(bnd.bcs[1])<:MatchedBC
        # (isempty(bnd.bcs[1]) & isempty(bnd.bcs[2])) ||
    # elseif typeof(bnd.bcs[1])<:MatchedBC
#
    # end

    i1 = latticeindex(lat,bnd.shape.start)
    i2 = latticeindex(lat,bnd.shape.stop)
    imin,imax = floor(Int,min(i1,i2)), floor(Int,max(i1,i2))
    x = lat[imin:imax]
    inds_keep = bnd.shape.(x)
    x = x[inds_keep]
    indices = CartesianIndex.(imin:imax)[inds_keep]
    ε = dielectric.(x)
    F = pump.(x)
    interior = trues(size(x))
    surface = falses(size(x))
    surface[1] = surface[end] = true
    bulk = interior .& .!surface
    corner = falses(size(x))
    nnm = (vcat(0,1:length(x)-1),)
    nnp = (vcat(2:length(x),0),)

    return Domain(type,name,bnd,lat,dielectric,pump,x,indices,ε,F,interior,bulk,surface,corner,nnm,nnp)
end
