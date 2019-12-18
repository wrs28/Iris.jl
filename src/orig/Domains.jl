# const SURFACE_BUFFER_FACTOR = 0/4 # abandoned for now

# """
#     Cavity(boundary,lattice,[dielectric,pump];align=false,[name]) -> domain
#     Cavity(lattice,boundary,[dielectric,pump];align=false,[name]) -> domain
# """
# Cavity(args...;kwargs...) = Domain(:Cavity,args...;kwargs...)
# """
#     Resonator(boundary,lattice,[dielectric,pump];align=false,[name]) -> domain
#     Resonator(lattice,boundary,[dielectric,pump];align=false,[name]) -> domain
# """
# Resonator(args...;kwargs...) = Domain(:Resonator,args...;kwargs...)
# """
#     Void(boundary,lattice;align=false,[name]) -> domain
#     Void(lattice,boundary;align=false,[name]) -> domain
# """
# Void(args...; align=false, kwargs...) = Domain(:Void,args...; align=align, kwargs...)
# """
#     Dielectric(lattice,boundary,[ε=1,F=0];align=false,[name]) -> domain
#     Dielectric(boundary,lattice,[ε=1,F=0];align=false,[name]) -> domain
# """
# Dielectric(args...;kwargs...) = Domain(:Dielectric,args...;kwargs...)
# Dielectric(bnd::Union{Boundary,Lattice},lat::Union{Boundary,Lattice},ε::Number,args...;kwargs...) = Domain(:Dielectric,bnd,lat,DielectricFunction(ε),args...;kwargs...)
# Dielectric(bnd::Union{Boundary,Lattice},lat::Union{Boundary,Lattice},ε,F::Number;kwargs...) = Domain(:Dielectric,bnd,lat,ε,PumpFunction(F);kwargs...)
# """
#     Resonator(boundary,lattice,[dielectric,pump];align=true,[name]) -> domain
#     Resonator(lattice,boundary,[dielectric,pump];align=true,[name]) -> domain
# """
# Waveguide(args...; kwargs...) = Domain(:Waveguide,args...; align=true, kwargs...)
# """
#     Lead(boundary,lattice,[dielectric,pump];align=true,[name]) -> domain
#     Lead(lattice,boundary,[dielectric,pump];align=true,[name]) -> domain
# """
# Lead(args...; kwargs...) = Domain(:Lead,args...; align=true, kwargs...)


# # polar spanning indices
# function generate_spanning_indices_polar(bnd::Boundary, lat::Lattice)
#     @assert typeof(bnd.shape)<:Union{Circle,Annulus} "polar coordinates only for circular outermost region"
#     @assert lat.x0==bnd.shape.x0 && lat.y0==bnd.shape.y0 "polar lattice origin $((lat.x0,lat.y0)) must be same as circle origin $((bnd.shape.x0,bnd.shape.y0))"
#     # i,j = get_lattice_index(lat,bnd.shape.x0,bnd.shape.y0)
#     # imin, imax = floor(Int,i), ceil(Int,i)
#     # jmin, jmax = floor(Int,j), ceil(Int,j)
#     # imin ≤ 0 ? imin = 0 : nothing
#     # jmin ≤ 0 ? jmin = 0 : nothing
#     imin = 0
#     imax = 1
#     jmin = 0
#     jmax = floor(Int,2π/lat.dθ)
#
#     flag_rmin = flag_rmax = flag_θmin = flag_θmax = true
#     init_flag = false
#     while flag_rmin || flag_rmax || flag_θmin|| flag_θmax
#         I = imin:imax
#         J = jmin:jmax
#         XYBool = BitArray(undef,2length(I)+2length(J))
#         count1 = 0
#         count2 = 0
#         temp = Array{Float64}(undef,2)
#         for i ∈ I
#             count2 += 1
#             X, Y = lat(i,J[1])
#             XYBool[1+count1+i-I[1]] = bnd.shape(X,Y)
#         end
#         (!init_flag || any(XYBool[(count1+1):count2])) && jmin≥1 ? (jmin -= 1; flag_θmin=true) : flag_θmin=false
#         !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
#         count1 = count2
#         for i ∈ I
#             count2 += 1
#             X,Y = lat(i,J[end])
#             XYBool[1+count1+i-I[1]] = bnd.shape(X,Y)
#         end
#         (!init_flag || any(XYBool[(count1+1):count2])) && (jmax+1)*lat.dθ≤2π ? (jmax += 1; flag_θmax=true) : flag_θmax=false
#         !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
#         count1 = count2
#         for j ∈ J
#             count2 += 1
#             X,Y = lat(I[1],j)
#             XYBool[1+count1+j-J[1]] = bnd.shape(X,Y)
#         end
#         (!init_flag || any(XYBool[(count1+1):count2])) && imin≥1 ? (imin -= 1; flag_rmin=true) : flag_rmin=false
#         !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
#         count1 = count2
#         for j ∈ J
#             count2 += 1
#             X,Y = lat(I[end],j)
#             XYBool[2length(I) + length(J) + j+1-J[1]] = bnd.shape(X,Y)
#         end
#         !init_flag || any(XYBool[(count1+1):count2]) ? (imax += 1; flag_rmax=true) : flag_rmax=false
#         !init_flag ? init_flag = any(XYBool[(count1+1):count2]) : nothing
#     end
#     return imin, imax, jmin, jmax
# end

# function generate_xy(bnd::Boundary,lat::Lattice,imin,imax,jmin,jmax)
#     I, J = imin:imax, jmin:jmax
#     X = Array{Float64}(undef,length(I),length(J))
#     Y = Array{Float64}(undef,length(I),length(J))
#     for k ∈ CartesianIndices(X)
#         i,j = k[1],k[2]
#         X[k],Y[k] = lat(I[i],J[j])
#     end
#     return X, Y
# end



# function generate_interior(bnd::Boundary,lat::Lattice,x::Array,y::Array)
#     interior = BitArray(undef,size(x)...)
#     for k ∈ CartesianIndices(x)
#         interior[k] = bnd.shape(x[k],y[k])
#     end
#     return interior
# end


# function generate_translations(lat::Lattice,imin,imax,jmin,jmax)
#     I, J = imin:imax, jmin:jmax
#     x1 = Array{Float64}(undef,length(I),length(J))
#     x2 = Array{Float64}(undef,length(I),length(J))
#     x3 = Array{Float64}(undef,length(I),length(J))
#     x4 = Array{Float64}(undef,length(I),length(J))
#     y1 = Array{Float64}(undef,length(I),length(J))
#     y2 = Array{Float64}(undef,length(I),length(J))
#     y3 = Array{Float64}(undef,length(I),length(J))
#     y4 = Array{Float64}(undef,length(I),length(J))
#     for k ∈ CartesianIndices(x1)
#         i,j = k[1],k[2]
#         x1[k],y1[k] = lat(I[i]-1,J[j])
#         x2[k],y2[k] = lat(I[i]+1,J[j])
#         x3[k],y3[k] = lat(I[i],J[j]-1)
#         x4[k],y4[k] = lat(I[i],J[j]+1)
#     end
#     return (x1,y1),(x2,y2),(x3,y3),(x4,y4)
# end


function generate_surface!(interior::BitArray,lat::Lattice,bnd::Boundary,t1,t2,t3,t4,x,y)
    surface = falses(size(interior)...)
    surface_temp = BitArray(undef,size(interior)...)
    for k ∈ CartesianIndices(interior)
        surface_temp[k] = interior[k] &&
            (!bnd.shape(t1[1][k],t1[2][k]) ||
             !bnd.shape(t2[1][k],t2[2][k]) ||
             !bnd.shape(t3[1][k],t3[2][k]) ||
             !bnd.shape(t4[1][k],t4[2][k])
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


function generate_neighbor_indices(interior,lat::Lattice)
    if lat.type==:Cartesian
        return generate_neighbor_indices_cartesian(interior)
    elseif lat.type==:Polar
        return generate_neighbor_indices_polar(interior)
    else
        throw(LatticeError(lat))
    end
end

function generate_neighbor_indices_cartesian(interior)
    NN = LinearIndices(interior) - reshape(cumsum(.!interior[:]),size(interior)...)
    nnxm = vcat(zeros(Int,1,size(interior,2)),NN[1:end-1,:])
    nnxp = vcat(NN[2:end,:],zeros(Int,1,size(interior,2)))
    nnym = hcat(zeros(Int,size(interior,1),1),NN[:,1:end-1])
    nnyp = hcat(NN[:,2:end],zeros(Int,size(interior,1),1))
    return nnxm,nnxp,nnym,nnyp
end

function generate_neighbor_indices_polar(interior)
    NN = LinearIndices(interior) - reshape(cumsum(.!interior[:]),size(interior)...)
    nnrm = vcat(zeros(Int,1,size(interior,2)),NN[1:end-1,:])
    nnrp = vcat(NN[2:end,:],zeros(Int,1,size(interior,2)))
    nnθm = hcat(NN[:,end],NN[:,1:end-1])
    nnθp = hcat(NN[:,2:end],NN[:,1])
    return nnrm,nnrp,nnθm,nnθp
end


function generate_corner(x,y,bnd::Boundary,lat::Lattice,surface,interior)
    if lat.type==:Cartesian
        return generate_corner_cartesian(x,y,bnd,lat,surface,interior)
    elseif lat.type==:Polar
        return generate_corner_polar(x,y,bnd,lat,surface,interior)
    else
        throw(LatticeError(lat))
    end
end


function generate_corner_cartesian(x,y,bnd::Boundary,lat::Lattice,surface,interior)
    corner = falses(size(surface)...)
    c = bnd.shape.corners
    dx = hypot(lat.dx,lat.dy)#*(1+SURFACE_BUFFER_FACTOR)
    if length(c)>0
        for k ∈ CartesianIndices(surface)
            if surface[k]
                d = map(z->hypot(x[k]-z[1],y[k]-z[2]),c)
                corner[k] = minimum(d) ≤ dx+eps(dx)
                X1, Y1 = x[k+CartesianIndex( 2, 0)], y[k+CartesianIndex( 2, 0)]
                X2, Y2 = x[k+CartesianIndex(-2, 0)], y[k+CartesianIndex(-2, 0)]
                X3, Y3 = x[k+CartesianIndex( 0, 2)], y[k+CartesianIndex( 0, 2)]
                X4, Y4 = x[k+CartesianIndex( 0,-2)], y[k+CartesianIndex( 0,-2)]
                count = sum(
                            (!bnd.shape(X1,Y1) && !interior[k+CartesianIndex( 1, 0)],
                             !bnd.shape(X2,Y2) && !interior[k+CartesianIndex(-1, 0)],
                             !bnd.shape(X3,Y3) && !interior[k+CartesianIndex( 0, 1)],
                             !bnd.shape(X4,Y4) && !interior[k+CartesianIndex( 0,-1)],)
                            )
                corner[k] = corner[k] && (count>1)
            end
        end
    end
    return corner
end

function generate_corner_polar(x,y,bnd,lat,surface,interior)
    return falses(size(surface)...)
end
