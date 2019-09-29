"""
    resonance_eigenproblem(simulation, k, ka=0, kb=0) -> A,B,σ

generate matrices from `simulation` to solve {∇²+εk²}ψ=0 in the form `A*v=λB*b` with shift `σ`
"""
function resonance_eigenproblem(sim::Simulation, k::Number, ka::Number=0, kb::Number=0)
    L = Laplacian(sim,k,ka,kb)
    return L.D, -L.M, k^2
end


"""
    resonance_eigenproblem(sim, k, ka=0, kb=0; η=0) -> A,B,σ

generate matrices from `simulation` to solve {∇²+εk²}u=-ηFk²u in the form `A*v=λB*b` with shift `σ
"""
function cf_eigenproblem(sim::Simulation, k::Number, ka::Number=0, kb::Number=0; η::Number=0)
    L = Laplacian(sim,k,ka,kb)
    return L.D+L.M*k^2, -L.Mf*k^2, η
end
