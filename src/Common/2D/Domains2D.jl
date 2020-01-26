module Domains2D


using ..Boundaries
using ..DielectricFunctions
using ..PumpFunctions
using ..Dispersions
using ..Lattices
using ..Points
using ..Shapes
using RecipesBase

import ..Symmetric, ..Unsymmetric
import LinearAlgebra.norm

import ..LatticeDomain

function LatticeDomain(
            boundary::Boundary{2},
            lattice::Lattice{2,Cartesian},
            n::Number = 1;
            type::Symbol = :generic,
            name::Symbol = :anonymous,
            fit::Bool = boundary.shape.ϕ ≈ lattice.ϕ)

    bnd = boundary
    lat = lattice

    if (typeof(bnd.shape) <: AbstractRectangle) && fit
        nx = ceil(Int,bnd.shape.a/lat.dx)
        ny = ceil(Int,bnd.shape.b/lat.dy)
        dx = bnd.shape.a/nx
        dy = bnd.shape.a/ny
        lat = Lattices.Cartesian(dx, dy; origin=bnd.shape.origin)
        CLASS = Symmetric()
    else
        CLASS = Unsymmetric()
    end

    return lattice_domain_classy(boundary,lattice,n,CLASS;type=type,name=name)
end


function LatticeDomain(
            boundary::Boundary{2},
            lattice::Lattice{2,Polar},
            n::Number = 1;
            type::Symbol = :generic,
            name::Symbol = :anonymous,
            fit::Bool = boundary.shape.origin ≈ lattice.origin)

    bnd = boundary
    lat = lattice

    if (typeof(bnd.shape) <: AbstractUniformDisk) && fit
        nr = ceil(Int,bnd.shape.r/lat.dr)
        nϕ = ceil(Int,2π/lat.dϕ)
        dr = bnd.shape.r/nr
        dϕ = 2π/nϕ
        lat = Lattices.Polar(dr, dϕ; origin=bnd.shape.origin)
        CLASS = Symmetric()
    else
        CLASS = Unsymmetric()
    end

    return lattice_domain_classy(boundary,lattice,n,CLASS;type=type,name=name)
end

function lattice_domain_classy(
            boundary::Boundary{2},
            lattice::Lattice{2},
            n::Number,
            ::CLASS;
            type::Symbol = :generic,
            name::Symbol = :anonymous
            ) where CLASS

    bnd = boundary
    lat = lattice

    ijminmax = _generate_spanning_indices(bnd,lat)
    p = _generate_points(bnd,lat,ijminmax...)
    interior = _generate_interior(bnd,lat,p)
    translations = _generate_translations(bnd,lat,ijminmax...)
    surface = _generate_surface!(interior,bnd,lat,p,translations...)
    nnxm, nnxp, nnym, nnyp = _generate_neighbor_indices(bnd,lat,interior)
    corner = _generate_corner(bnd,lat,p,surface,interior)

    inds_keep = findall(interior)
    x = p[inds_keep]
    interior = interior[inds_keep]
    surface = surface[inds_keep]
    bulk = interior .& .!surface
    corner = corner[inds_keep]
    nnxm = nnxm[inds_keep]
    nnxp = nnxp[inds_keep]
    nnym = nnym[inds_keep]
    nnyp = nnyp[inds_keep]
    indices = map(z->z+CartesianIndex(ijminmax[1]-1,ijminmax[3]-1),inds_keep)

    return LatticeDomain{CLASS}(bnd, lat, n, x, indices, interior,
                            bulk, surface, corner, (nnxm,nnym), (nnxp,nnyp); type=type, name=name)
end


function Base.propertynames(::LatticeDomain{2}, private=false)
    if private
        return fieldnames(LatticeDomain)
    else
        return (:boundary, :lattice, :shape, :n, :ε, :type, :name, :x, :interior, :surface, :corner, :bulk)
    end
end

################################################################################
# DOMAIN BUILDING UTITLITIES

# cartesian spanning indices
function _generate_spanning_indices(bnd::Boundary{2}, lat::Lattice{2,Cartesian})
    i,j = latticeindex(lat,bnd.shape.origin)
    imin, imax = floor(Int,i), ceil(Int,i)
    jmin, jmax = floor(Int,j), ceil(Int,j)

    flag_xmin = flag_xmax = flag_ymin = flag_ymax = true
    init_flag = false
    while flag_ymin || flag_ymax || flag_xmin|| flag_xmax
        I = imin:imax
        J = jmin:jmax
        PBool = BitArray(undef,2length(I)+2length(J))
        count1 = 0
        count2 = 0
        temp = Vector{Float64}(undef,2)
        for i ∈ I
            count2 += 1
            P = lat[i,J[1]]
            PBool[1+count1+i-I[1]] = bnd.shape(P)
        end
        !init_flag || any(PBool[(count1+1):count2]) ? (jmin -= 1; flag_ymin=true) : flag_ymin=false
        !init_flag ? init_flag = any(PBool[(count1+1):count2]) : nothing
        count1 = count2
        for i ∈ I
            count2 += 1
            P = lat[i,J[end]]
            PBool[1+count1+i-I[1]] = bnd.shape(P)
        end
        !init_flag || any(PBool[(count1+1):count2]) ? (jmax += 1; flag_ymax=true) : flag_ymax=false
        !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
        count1 = count2
        for j ∈ J
            count2 += 1
            P = lat[I[1],j]
            PBool[1+count1+j-J[1]] = bnd.shape(P)
        end
        !init_flag || any(PBool[(count1+1):count2]) ? (imin -= 1; flag_xmin=true) : flag_xmin=false
        !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
        count1 = count2
        for j ∈ J
            count2 += 1
            P = lat[I[end],j]
            PBool[2length(I) + length(J) + j+1-J[1]] = bnd.shape(P)
        end
        !init_flag || any(PBool[(count1+1):count2]) ? (imax += 1; flag_xmax=true) : flag_xmax=false
        !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
    end
    return imin-2, imax+2, jmin-2, jmax+2
end

# polar spanning indices
function _generate_spanning_indices(bnd::Boundary, lat::Lattice{2,Polar})
    @assert typeof(bnd.shape)<:Union{Circle,Annulus} "polar coordinates only for circular outermost region"
    @assert lat.x0==bnd.shape.x0 && lat.y0==bnd.shape.y0 "polar lattice origin $((lat.x0,lat.y0)) must be same as circle origin $((bnd.shape.x0,bnd.shape.y0))"
    # i,j = get_lattice_index(lat,bnd.shape.x0,bnd.shape.y0)
    # imin, imax = floor(Int,i), ceil(Int,i)
    # jmin, jmax = floor(Int,j), ceil(Int,j)
    # imin ≤ 0 ? imin = 0 : nothing
    # jmin ≤ 0 ? jmin = 0 : nothing
    imin = 0
    imax = 1
    jmin = 0
    jmax = floor(Int,2π/lat.dθ)

    flag_rmin = flag_rmax = flag_θmin = flag_θmax = true
    init_flag = false
    while flag_rmin || flag_rmax || flag_θmin|| flag_θmax
        I = imin:imax
        J = jmin:jmax
        XYBool = BitArray(undef,2length(I)+2length(J))
        count1 = 0
        count2 = 0
        temp = Array{Float64}(undef,2)
        for i ∈ I
            count2 += 1
            X, Y = lat(i,J[1])
            XYBool[1+count1+i-I[1]] = bnd.shape(X,Y)
        end
        (!init_flag || any(XYBool[(count1+1):count2])) && jmin≥1 ? (jmin -= 1; flag_θmin=true) : flag_θmin=false
        !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
        count1 = count2
        for i ∈ I
            count2 += 1
            X,Y = lat(i,J[end])
            XYBool[1+count1+i-I[1]] = bnd.shape(X,Y)
        end
        (!init_flag || any(XYBool[(count1+1):count2])) && (jmax+1)*lat.dθ≤2π ? (jmax += 1; flag_θmax=true) : flag_θmax=false
        !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
        count1 = count2
        for j ∈ J
            count2 += 1
            X,Y = lat(I[1],j)
            XYBool[1+count1+j-J[1]] = bnd.shape(X,Y)
        end
        (!init_flag || any(XYBool[(count1+1):count2])) && imin≥1 ? (imin -= 1; flag_rmin=true) : flag_rmin=false
        !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
        count1 = count2
        for j ∈ J
            count2 += 1
            X,Y = lat(I[end],j)
            XYBool[2length(I) + length(J) + j+1-J[1]] = bnd.shape(X,Y)
        end
        !init_flag || any(XYBool[(count1+1):count2]) ? (imax += 1; flag_rmax=true) : flag_rmax=false
        !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
    end
    return imin, imax, jmin, jmax
end

function _generate_points(bnd::Boundary{2},lat::Lattice{2},imin,imax,jmin,jmax)
    I, J = imin:imax, jmin:jmax
    return lat[[i for i ∈ I, j ∈ J],[j for i ∈ I, j ∈ J]]
end

function _generate_interior(bnd::Boundary{2},lat::Lattice{2},p)
    interior = BitArray(undef,size(p)...)
    for k ∈ CartesianIndices(p) interior[k] = bnd.shape(p[k]) end
    return interior
end

function _generate_translations(bnd::Boundary{2},lat::Lattice{2},imin,imax,jmin,jmax)
    I, J = imin:imax, jmin:jmax
    p1 = Matrix{Point{2,Cartesian}}(undef,length(I),length(J))
    p2 = Matrix{Point{2,Cartesian}}(undef,length(I),length(J))
    p3 = Matrix{Point{2,Cartesian}}(undef,length(I),length(J))
    p4 = Matrix{Point{2,Cartesian}}(undef,length(I),length(J))
    for k ∈ CartesianIndices(p1)
        i,j = k[1],k[2]
        p1[k] = Cartesian(lat[I[i]-1,J[j]])
        p2[k] = Cartesian(lat[I[i]+1,J[j]])
        p3[k] = Cartesian(lat[I[i],J[j]-1])
        p4[k] = Cartesian(lat[I[i],J[j]+1])
    end
    return p1, p2, p3, p4
end

function _generate_surface!(interior::BitArray,bnd::Boundary{2},lat::Lattice{2},p,t1,t2,t3,t4)
    surface = falses(size(interior)...)
    surface_temp = BitArray(undef,size(interior)...)
    for k ∈ CartesianIndices(interior)
        surface_temp[k] =   interior[k] &&
                            (!bnd.shape(t1[k]) ||
                             !bnd.shape(t2[k]) ||
                             !bnd.shape(t3[k]) ||
                             !bnd.shape(t4[k])
                            )
        if surface_temp[k]
            # n,t,d = normal_distance(bnd.shape,x[k],y[k])
            # if d[1]<min(lat.dx*SURFACE_BUFFER_FACTOR,lat.dy*SURFACE_BUFFER_FACTOR)
                # interior[k] = false
                # surface[k] = false
                # surface[k+CartesianIndex( 1, 0)]=true
                # surface[k+CartesianIndex(-1, 0)]=true
                # surface[k+CartesianIndex( 0, 1)]=true
                # surface[k+CartesianIndex( 0,-1)]=true
            # else
                surface[k] = true
            # end
        end
    end
    return surface
end

function _generate_neighbor_indices(bnd::Boundary{2},lat::Lattice{2,Cartesian},interior)
    NN = LinearIndices(interior) - reshape(cumsum(.!interior[:]),size(interior)...)
    nnxm = vcat(zeros(Int,1,size(interior,2)),NN[1:end-1,:])
    nnxp = vcat(NN[2:end,:],zeros(Int,1,size(interior,2)))
    nnym = hcat(zeros(Int,size(interior,1),1),NN[:,1:end-1])
    nnyp = hcat(NN[:,2:end],zeros(Int,size(interior,1),1))
    return nnxm,nnxp,nnym,nnyp
end

function _generate_neighbor_indices(bnd::Boundary{2},lat::Lattice{2,Polar},interior)
    NN = LinearIndices(interior) - reshape(cumsum(.!interior[:]),size(interior)...)
    nnrm = vcat(zeros(Int,1,size(interior,2)),NN[1:end-1,:])
    nnrp = vcat(NN[2:end,:],zeros(Int,1,size(interior,2)))
    nnθm = hcat(NN[:,end],NN[:,1:end-1])
    nnθp = hcat(NN[:,2:end],NN[:,1])
    return nnrm,nnrp,nnθm,nnθp
end

function _generate_corner(bnd::Boundary{2},lat::Lattice{2,Cartesian},p,surface,interior)
    corner = falses(size(surface)...)
    c = bnd.shape.corners
    dx = hypot(lat.dx,lat.dy)#*(1+SURFACE_BUFFER_FACTOR)
    if length(c)>0
        for k ∈ CartesianIndices(surface)
            if surface[k]
                d = map(z->norm(p[k]-z),c)
                corner[k] = minimum(d) ≤ dx+eps(dx)
                P1 = p[k+CartesianIndex( 2, 0)]
                P2 = p[k+CartesianIndex(-2, 0)]
                P3 = p[k+CartesianIndex( 0, 2)]
                P4 = p[k+CartesianIndex( 0,-2)]
                count = sum(
                            (!bnd.shape(P1) && !interior[k+CartesianIndex( 1, 0)],
                             !bnd.shape(P2) && !interior[k+CartesianIndex(-1, 0)],
                             !bnd.shape(P3) && !interior[k+CartesianIndex( 0, 1)],
                             !bnd.shape(P4) && !interior[k+CartesianIndex( 0,-1)],)
                            )
                corner[k] = corner[k] && (count>1)
            end
        end
    end
    return corner
end

_generate_corner(bnd::Boundary{2},lat::Lattice{2,Polar},p,surface,interior) = falses(size(surface)...)

################################################################################
# PLOTTING

import ...Common.MARKERSIZE_SCALE
import ...Common.MARKERSHAPE

@recipe function f(d::LatticeDomain{2})
    aspect_ratio --> 1
    legend --> false
    @series begin
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(d.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
        d.x[d.bulk]
    end
    @series begin
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(d.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
        d.x[d.surface]
    end
    @series begin
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(d.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
        d.x[d.corner]
    end
    # @series d.boundary
end

end

using .Domains2D
