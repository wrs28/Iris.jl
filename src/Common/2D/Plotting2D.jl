import ..Common.MARKERSIZE_SCALE
import ..Common.MARKERSHAPE

# Simulation
@recipe function f(sim::Simulation{2})
    n = sqrt.(sim.ε[:])

    layout := (1,3)
    legend --> false

    @series begin
		title --> L"{\rm Real}(n)"
        subplot := 1
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(sim.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
		seriescolor --> :sequential
		aspectratio --> 1
		zcolor --> real(n)
		clims --> (1,3)
        sim.x
    end
	@series begin
		title --> L"{\rm Imag}(n)"
        subplot := 2
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(sim.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
		seriescolor --> :diverging
		aspectratio --> 1
		zcolor --> imag(n)
		clims --> (-1,1)
        sim.x
    end
	@series begin
		title --> L"{\rm Pump\ }F"
        subplot := 3
        seriestype --> :scatter
        markersize --> MARKERSIZE_SCALE/sqrt(length(sim.x))
        markerstrokealpha --> 0
        shape --> MARKERSHAPE
		seriescolor --> :diverging
		aspectratio --> 1
		subplot := 3
		zcolor --> sim.F[:]
		clims --> (-1,1)
        sim.x
    end
end
#
# # Electric Field
# @recipe function f(e::ElectricField{1}, by::Function)
#     x = map(a->a.vec[1],e.pos)
#     perm = sortperm(x)
# 	x = x[perm]
#
# 	upper_lim = maximum(abs∘by,e.val)
# 	if by ∈ (abs,abs2)
# 		lower_lim = zero(upper_lim)
# 	else
# 		lower_lim = -upper_lim
# 	end
#
# 	m = size(e,2)
#     layout --> (3,m)
#     legend --> false
# 	ylims --> (lower_lim,upper_lim)
# 	for μ ∈ 1:m
# 	    @series begin
# 			title --> latexstring("{\\rm Field}\\ $μ")
# 	        ylabel --> L"E_x"
# 	        subplot := 0m+μ
# 	        x, by.(e(μ).x[perm])
# 	    end
# 	    @series begin
# 	        ylabel --> L"E_y"
# 	        subplot := 1m+μ
# 	        x, by.(e(μ).y[perm])
# 	    end
# 	    @series begin
# 			xlabel --> L"x"
# 	        ylabel --> L"E_z"
# 	        subplot := 2m+μ
# 	        x, by.(e(μ).z[perm])
# 	    end
# 	end
# end
#
#
# @recipe function f(sim::Simulation{1},e::ElectricField{1}, by::Function)
#     x = map(a->a.vec[1],e.pos)
#     perm = sortperm(x)
#     x = x[perm]
#
# 	legend --> false
#
# 	m = 1+size(e,2)
# 	layout --> (3,m)
#
# 	n = sqrt.(sim.ε[perm])
# 	@series begin
# 		title --> L"{\rm Scattering\ Structure}"
# 		ylabel --> L"{\rm real}(n)"
#         subplot := 0m+1
#         x, real(n)
#     end
#     @series begin
# 		ylabel --> L"{\rm imag}(n)"
#         subplot := 1m+1
#         x, imag(n)
#     end
#     @series begin
# 		xlabel --> L"x"
# 		ylabel --> L"{\rm pump\ F}"
# 		subplot := 2m+1
#         x, sim.F[perm]
#     end
#
# 	upper_lim = maximum(abs∘by,e.val)
# 	if by ∈ (abs,abs2)
# 		lower_lim = zero(upper_lim)
# 	else
# 		lower_lim = -upper_lim
# 	end
# 	ylims --> (lower_lim,upper_lim)
#
# 	for μ ∈ 1:m-1
# 	    @series begin
# 			title --> latexstring("{\\rm Field}\\ $μ")
# 	        ylabel --> L"E_x"
# 	        subplot := 0m+μ+1
# 	        x, by.(e(μ).x[perm])
# 	    end
# 	    @series begin
# 	        ylabel --> L"E_y"
# 	        subplot := 1m+μ+1
# 	        x, by.(e(μ).y[perm])
# 	    end
# 	    @series begin
# 			xlabel --> L"x"
# 	        ylabel --> L"E_z"
# 	        subplot := 2m+μ+1
# 	        x, by.(e(μ).z[perm])
# 	    end
#
# 	end
# end
