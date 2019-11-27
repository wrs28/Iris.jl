using IterativeSolvers

fnames = (:cg,:minres,:gmres,:idrs,:bicgstabl,:chebyshev,:jacobi,:gauss_seidel,:sor,:ssor,:lsmr,:lsqr)
for fname ∈ fnames
    @eval begin
        function IterativeSolvers.$(fname)(mls::MaxwellLS, args...; kwargs...)
            $(Symbol(fname,"!"))(view(mls.solution.tot.values,:,1), mls.maxwell.A, mls.j, args...; kwargs...)
            return mls.solution
        end
        function IterativeSolvers.$(Symbol(fname,"!"))(mls::MaxwellLS, args...; kwargs...)
            _, hist = $(Symbol(fname,"!"))(view(mls.solution.tot.values,:,1), mls.maxwell.A, reshape(Array(mls.j),length(mls.j),:), args...; kwargs..., log=true)
            for i ∈ eachindex(mls.equivalent_source.field)
                k = mod1(i,length(mls.j)÷3)
                mls.solution.incident[i] = mls.equivalent_source.field[i]*mls.equivalent_source.mask(mls.equivalent_source.field.pos[k])
                mls.solution.scattered[i] = mls.solution.total[i] - mls.solution.incident[i]
            end
            mls.solved[] = true
            mls.converged[] = hist.isconverged
            return nothing
        end
    end
end
