_bls_by_shape(bls,s::Square) = fix_bls_by_shape!(bls,Rectangle(s))
function _bls_by_shape(bls,shape::Rectangle)
	stop = shape.x0-shape.a/2
	start = stop + bls[1].depth
	bls[1] = bls[1](startx=start, stopx=stop)

	stop = shape.x0+shape.a/2
	start = stop - bls[2].depth
	bls[2] = bls[2](startx=start, stopx=stop)

	stop = shape.y0-shape.b/2
	start = stop + bls[3].depth
	bls[3] = bls[3](starty=start, stopy=stop)

	stop = shape.y0+shape.b/2
	start = stop - bls[4].depth
	bls[4] = bls[4](starty=start, stopy=stop)
end
function _bls_by_shape(bls,c::Circle)
	stop=c.R
	start=stop-bls[1].depth
	bls[1] = bls[1](startx=start,stopx=stop)
end



#
# """
# 	plot(::Boundary)
#
# PML's plotted in blue patches (b/c absorbing, so cold), cPML's in red (b/c emitting, so hot)
#
# DirichletBC are solid black lines, NeumannBC dashed, FloquetBC dashdot, MatchedBC/noBC nothing
# """
# @recipe function f(bnd::Boundary)
# 	aspect_ratio --> 1
# 	legend --> false
#
# 	@series begin
# 		fillcolor --> BOUNDARY_COLOR
# 		color --> BOUNDARY_COLOR
# 		markershape := :none
# 		bnd.shape
# 	end
# 	@series (bnd,1)
# 	@series (bnd,1,1)
# end
# # plot boundary layers
# @recipe function f(bnd::Boundary,_1::Int)
# 	shape = bnd.shape
# 	cosθ = shape.cosθ
# 	sinθ = shape.sinθ
# 	typeof(shape)<:Square ? shape = Rectangle(shape) : nothing
# 	for i ∈ eachindex(bnd.bls)
# 		bls = bnd.bls[i]
# 		depth = bnd.bls[i].depth
# 		if typeof(shape)<:Rectangle
# 			if i==1
# 				x = -(shape.a-depth)/2
# 				y = 0
# 				bl_shape = Rectangle(depth,shape.b,shape.x0+cosθ*x-sinθ*y,shape.y0+sinθ*x+cosθ*y,shape.θ)
# 			elseif i==4
# 				x = 0
# 				y = (shape.b-depth)/2
# 				bl_shape = Rectangle(shape.a,depth,shape.x0+cosθ*x-sinθ*y,shape.y0+sinθ*x+cosθ*y,shape.θ)
# 			elseif i==2
# 				x = +(shape.a-depth)/2
# 				y = 0
# 				bl_shape = Rectangle(depth,shape.b,shape.x0+cosθ*x-sinθ*y,shape.y0+sinθ*x+cosθ*y,shape.θ)
# 			else
# 				x = 0
# 				y = -(shape.b-depth)/2
# 				bl_shape = Rectangle(shape.a,depth,shape.x0+cosθ*x-sinθ*y,shape.y0+sinθ*x+cosθ*y,shape.θ)
# 			end
# 		elseif typeof(shape)<:Circle
# 			bls = bnd.bls[1]
# 			depth = bls.depth
# 			bl_shape = Annulus(shape.R-depth,shape.R,shape.x0,shape.y0)
# 		end
# 		if typeof(bls)<:AbstractComplexBL
# 			@series begin
# 				alpha --> 0
# 				typeof(bls)<:PML ? fillcolor --> BOUNDARY_PML_COLOR : fillcolor --> BOUNDARY_CPML_COLOR
# 				typeof(bls)<:PML ? color --> BOUNDARY_PML_COLOR : color --> BOUNDARY_CPML_COLOR
# 				bl_shape
# 			end
# 		end
# 	end
# end
# # plot boundary conditions
# @recipe function f(bnd::Boundary,_1::Int,_2::Int)
# 	shape = bnd.shape
# 	for i ∈ eachindex(bnd.bcs)
# 		bcs = bnd.bcs[i]
# 		if typeof(shape)<:Rectangle
# 			if i==1
# 				x = [-1,-1]*shape.a/2
# 				y = [-1,+1]*shape.b/2
# 			elseif i==2
# 				x = [-1,+1]*shape.a/2
# 				y = [+1,+1]*shape.b/2
# 			elseif i==3
# 				x = [+1,+1]*shape.a/2
# 				y = [-1,+1]*shape.b/2
# 			else
# 				x = [-1,+1]*shape.a/2
# 				y = [-1,-1]*shape.b/2
# 			end
# 			if typeof(bcs)<:AbstractLocalBC
# 				@series begin
# 					if typeof(bcs)<:DirichletBC
# 						ls --> BOUNDARY_DIRICHLET_LINETYPE
# 					elseif typeof(bcs)<:NeumannBC
# 						ls --> BOUNDARY_NEUMANN_LINETYPE
# 					else
# 						ls --> BOUNDARY_FLOQUET_MATCHED_LINETYPE
# 					end
# 					seriestype --> :line
# 					markershape := :none
# 					color --> BOUNDARY_BC_COLOR
# 					shape.x0 .+ shape.cosθ*x-shape.sinθ*y, shape.y0 .+ shape.sinθ*x+shape.cosθ*y
# 				end
# 			end
# 		else
# 			if typeof(bcs)<:Union{AbstractLocalBC,FloquetBC}
# 				@series begin
# 					if typeof(bcs)<:DirichletBC
# 						ls --> BOUNDARY_DIRICHLET_LINETYPE
# 					elseif typeof(bcs)<:NeumannBC
# 						ls --> BOUNDARY_NEUMANN_LINETYPE
# 					else
# 						ls --> BOUNDARY_FLOQUET_MATCHED_LINETYPE
# 					end
# 					color --> BOUNDARY_BC_COLOR
# 					markershape := :none
# 					fillalpha --> 0
# 					alpha --> 1
# 					shape
# 				end
# 			end
# 		end
# 	end
# end
#
