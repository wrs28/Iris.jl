function _getdimensionalproperty(s::ScalarField{2}, sym::Symbol)
    n,m = size(getfield(s,:values))
    if sym==:left
        start_ind = getfield(s,:start_inds)[1]
        values = getfield(s,:values)
        inds = [start_ind]
        p = values[inds,:]
    elseif sym==:right
        stop_ind = getfield(s,:stop_inds)[1]
        values = getfield(s,:values)
        inds = [stop_ind]
        p = values[inds,:]
    elseif sym==:bottom
        start_ind = getfield(s,:start_inds)[2]
        values = getfield(s,:values)
        inds = [start_ind]
        p = values[inds,:]
    elseif sym==:top
        stop_ind = getfield(s,:stop_inds)[2]
        values = getfield(s,:values)
        inds = [stop_ind]
        p = values[inds,:]
    else
        return nothing
    end
end

function _getdimensionalproperty(e::ElectricField{2}, sym::Symbol)
    n,m = size(getfield(e,:values))
    if sym==:left
        start_ind = getfield(e,:start_inds)[1]
        values = getfield(e,:values)
        inds = [0n÷3+start_ind,1n÷3+start_ind,2n÷3+start_ind]
        p = values[inds,:]
    elseif sym==:right
        stop_ind = getfield(e,:stop_inds)[1]
        values = getfield(e,:values)
        inds = [0n÷3+stop_ind,1n÷3+stop_ind,2n÷3+stop_ind]
        p = values[inds,:]
    elseif sym==:bottom
        start_ind = getfield(e,:start_inds)[2]
        values = getfield(e,:values)
        inds = [0n÷3+start_ind,1n÷3+start_ind,2n÷3+start_ind]
        p = values[inds,:]
    elseif sym==:top
        stop_ind = getfield(e,:stop_inds)[2]
        values = getfield(e,:values)
        inds = [0n÷3+stop_ind,1n÷3+stop_ind,2n÷3+stop_ind]
        p = values[inds,:]
    else
        return nothing
    end
end

function Base.propertynames(::ScalarField{2}, private=false)
    if private
        return fieldnames(ScalarField)
    else
        return (:positions,:values,:left,:right,:bottom,:top)
    end
end

function Base.propertynames(::ElectricField{2}, private=false)
    if private
        return fieldnames(ElectricField)
    else
        return (:positions,:values,:x,:y,:z,:left,:right,:bottom,:top)
    end
end


# """
#     interpolate(::ElectricField{1}, interpmode, gridstyle)
#
# See documentation for `interpolate` from [`Interpolations.jl`](https://github.com/JuliaMath/Interpolations.jl)
# """
# Interpolations.interpolate(e::ElectricField{1}, interpmode, gridstyle) =
#     map(i->map(s->interpolate(getproperty(e(i),s), interpmode, gridstyle),(:x,:y,:z)),1:size(e,2))
#
# fnames = (:LinearInterpolation,:CubicSplineInterpolation)
# for fn ∈ fnames
#     @eval begin
#         """
#             $($fn)(::ElectricField{1}; kwargs...) -> interp
#
#         See documentation for `$($fn)` from [`Interpolations.jl`](https://github.com/JuliaMath/Interpolations.jl)
#         """
#         function Interpolations.$fn(e::ElectricField{1}; kwargs...)
#             x = map(p->p.x,e.pos)
#             perm = sortperm(x)
#             xs = LinRange(e.starts[1],e.stops[1],length(e.pos))
#             itps = map(i->InterpolatedElectricField(e.pos,map(s->$fn(xs,getproperty(e(i),s)[perm]; kwargs...), (:x,:y,:z))...), ntuple(identity,size(e,2)))
#             if length(itps)==1
#                 return itps[1]
#             else
#                 return itps
#             end
#         end
#     end
# end
