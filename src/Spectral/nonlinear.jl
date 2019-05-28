"""
    resonance_nonlinear_eigenproblem(sim, k, ka=0, kb=0; kwargs...)
"""
function resonance_nonlinear_eigenproblem(sim::Simulation, k::Number, ka::Number=0, kb::Number=0; kwargs...)
    N,brgs = sim.dis.N,(sim.dis.coordinate_system,sim.bnd.∂Ω,sim.bnd.bc,sim.bnd.bl)
    (∇²,fs,s),defs = laplacian(N,brgs...,k,ka,kb,0)
    for i ∈ eachindex(∇²)
        ∇²[i] = (∇²[i] + SparseMatrixCSC(transpose(∇²[i])))/2
    end
    B = s*spdiagm(0=>sim.sys.ε[:])
    push!(∇²,B)
    push!(fs,k->k^2)
    return SPMF_NEP(∇²,fs; check_consistency=false, kwargs...)
end
