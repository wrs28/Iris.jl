# function green_fn(sim::Simulation{1},ω::Number,x::Real)
#     inds = IrosBase.Tessellations.get_surrounding_sites(sim.tessellation,x,y)
#     weights = IrosBase.Simulations.generate_weights(x,y,sim.x[inds],sim.y[inds])
#     j = zeros(sim.n)
#     x = sim.x[inds]; y = sim.y[inds]
#     v1 = [x[2]-x[1],y[2]-y[1]]
#     v2 = [x[3]-x[2],y[3]-y[2]]
#     cosθ = (norm(v1+v2)^2 - norm(v1)^2 - norm(v2)^2)/2/norm(v1)/norm(v2)
#     sinθ = sqrt(1-cosθ^2)
#     area = norm(v1)*norm(v2)*sinθ/2
#     j[inds] = weights./area
#     return solve_scatter(sim,k,j)
# end

function scattering(
            sim::Simulation{1},
            ω::Number,
            j::AbstractArray;
            ky::Number=0,
            kz::Number=0)

    A,B = maxwell_lep(sim;ky=ky,kz=kz)(ω)
    return scattering_core(A-B*ω^2,j)
end
#
#
# function solve_scatter(sim::Simulation,kx::Number,ky::Number,shape::AbstractShape,f::Function)
#     k = hypot(kx,ky)
#     x = sim.x[sim.interior]
#     y = sim.y[sim.interior]
#     L = Laplacian(sim,k,0,0)
#     φ = (f.(kx,ky,x,y).*shape.(x,y))
#     j = (L.D+sparse(I*k^2,sim.n,sim.n))*φ
#     A = L.D + L.M*k^2
#     scattered_field = solve_scatter_core(A,-j)
#     total_field = scattered_field + φ
#     return total_field, scattered_field, φ
# end
# solve_scatter(sim::Simulation,k::Vector,shape::AbstractShape,f::Function) = solve_scatter(sim,k[1],k[2],shape,f)
