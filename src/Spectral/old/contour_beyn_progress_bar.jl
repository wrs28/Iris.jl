NonlinearEigenproblems.contour_beyn(nep::NEP,display::Bool;params...)=contour_beyn(ComplexF64,nep,display;params...)

# the same as NonlinearEigenproblems.contour_beyn but compatible with ProgressMeter, notice extra display argument to distinguish from main method
function NonlinearEigenproblems.contour_beyn(::Type{T},
            nep::NEP,
            display::Bool;
            tol::Real=sqrt(eps(real(T))), # Note tol is quite high for this method
            σ::Number=zero(complex(T)),
            displaylevel::Integer=0,
            linsolvercreator::Function=backslash_linsolvercreator,
            neigs::Integer=2, # Number of wanted eigvals
            k::Integer=neigs+1, # Columns in matrix to integrate
            radius::Union{Real,Tuple,Array}=1, # integration radius
            quad_method::Symbol=:ptrapz, # which method to run. :quadg, :quadg_parallel, :quadgk, :ptrapz
            N::Integer=1000,  # Nof quadrature nodes
            errmeasure::Function = NonlinearEigenproblems.default_errmeasure(nep::NEP),
            sanity_check=true,
            rank_drop_tol=tol # Used in sanity checking
            )where{T<:Number}

    # Geometry
    length(radius)==1 ? radius=(radius,radius) : nothing
    g(t) = complex(radius[1]*cos(t),radius[2]*sin(t)) # ellipse
    gp(t) = complex(-radius[1]*sin(t),radius[2]*cos(t)) # derivative

    n=size(nep,1);


    if (k>n)
        error("Cannot compute more eigenvalues than the size of the NEP with contour_beyn() k=",k," n=",n);

    end
    if (k<=0)
        error("k must be positive, k=",k,
              neigs==typemax(Int) ? ". The kwarg k must be set if you use neigs=typemax" : ".")
    end


    Random.seed!(10); # Reproducability
    Vh=Array{T,2}(randn(real(T),n,k)) # randn only works for real


    function local_linsolve(λ::TT,V::Matrix{TT}) where {TT<:Number}
#        @ifd(print("."))
        local M0inv::LinSolver = linsolvercreator(nep,λ+σ);
        # This requires that lin_solve can handle rectangular
        # matrices as the RHS
        return lin_solve(M0inv,V);
    end

    # Constructing integrands
    Tv(λ) = local_linsolve(T(λ),Vh)
    f(t) = Tv(g(t))*gp(t)
    # @ifd(print("Computing integrals"))

    pg = Progress(N; dt=.1, desc="Contour integration...")
    local A0,A1
    if (quad_method == :quadg_parallel)
        # @ifd(print(" using quadg_parallel"))
        error("disabled");
    elseif (quad_method == :quadg)
        # @ifd(print(" using quadg"))
        error("disabled");
    elseif (quad_method == :ptrapz)
        # @ifd(print(" using ptrapz"))
        (A0,A1)=ptrapz(f,g,0,2*pi,N,pg);
    elseif (quad_method == :ptrapz_parallel)
        # @ifd(print(" using ptrapz_parallel"))
        channel = RemoteChannel(()->Channel{Bool}(pg.n+1),1)
        @sync begin
            @async begin
                (A0,A1)=ptrapz_parallel(f,g,0,2*pi,N,channel);
                put!(channel,false)
            end
            @async while take!(channel)
                next!(pg)
            end
        end
    elseif (quad_method == :quadgk)
        error("disabled");
    else
        error("Unknown quadrature method:"*String(quad_method));
    end

    # Don't forget scaling
    A0[:,:] = A0 ./(2im*pi);
    A1[:,:] = A1 ./(2im*pi);

    # @ifd(print("Computing SVD prepare for eigenvalue extraction "))
    V,S,W = svd(A0)
    p = count( S/S[1] .> rank_drop_tol);
    # @ifd(println(" p=",p));

    V0 = V[:,1:p]
    W0 = W[:,1:p]
    B = (copy(V0')*A1*W0) * Diagonal(1 ./ S[1:p])

    # Extract eigenval and eigvec approximations according to
    # step 6 on page 3849 in the reference
    # @ifd(println("Computing eigenvalues "))
    λ,VB=eigen(B)
    λ[:] = λ .+ σ
    # @ifd(println("Computing eigenvectors "))
    V = V0 * VB;
    for i = 1:p
        normalize!(V[:,i]);
    end

    if (!sanity_check)
        sorted_index = sortperm(map(x->abs(σ-x), λ));
        inside_bool = (real(λ[sorted_index].-σ)/radius[1]).^2 + (imag(λ[sorted_index].-σ)/radius[2]).^2 .≤ 1
        inside_perm = sortperm(.!inside_bool)
        return (λ[sorted_index[inside_perm]],V[:,sorted_index[inside_perm]])
    end

    # Compute all the errors
    errmeasures=zeros(real(T),p);
    # for i = 1:p
        # errmeasures[i]=errmeasure(λ[i],V[:,i]);
    # end

    good_index=findall(errmeasures .< tol);

    # Index vector for sorted to distance from σ
    sorted_good_index=
       good_index[sortperm(map(x->abs(σ-x), λ[good_index]))];

    # check that eigenvalues are inside contour, move them to the end if they are outside
   inside_bool = (real(λ[sorted_good_index].-σ)/radius[1]).^2 + (imag(λ[sorted_good_index].-σ)/radius[2]).^2 .≤ 1
   if any(.!inside_bool)
       if neigs ≤ sum(inside_bool)
           @warn "found $(sum(.!inside_bool)) evals outside contour, $p inside. all $neigs returned evals inside, but possibly inaccurate. try increasing N, decreasing tol, or changing radius"
       else
           @warn "found $(sum(.!inside_bool)) evals outside contour, $p inside. last $(neigs-sum(inside_bool)) returned evals outside contour. possible inaccuracy. try increasing N, decreasing tol, or changing radius"
       end
   end
   sorted_good_inside_perm = sortperm(.!inside_bool)

    # Remove all eigpairs not sufficiently accurate
    # and potentially remove eigenvalues if more than neigs.
    local Vgood,λgood
    if( size(sorted_good_index,1) > neigs)
        # @ifd(println("Removing unwanted eigvals: neigs=",neigs,"<",size(sorted_good_index,1),"=found_eigvals"))
        Vgood=V[:,sorted_good_index[sorted_good_inside_perm][1:neigs]];
        λgood=λ[sorted_good_index[sorted_good_inside_perm][1:neigs]];
    else
        Vgood=V[:,sorted_good_index[sorted_good_inside_perm]];
        λgood=λ[sorted_good_index[sorted_good_inside_perm]];
    end

    if (p==k)
       @warn "Rank-drop not detected, your eigvals may be correct, but the algorithm cannot verify. Try to increase k. This warning can be disabled with `sanity_check=false`." S
    end

    if (size(λgood,1)<neigs  && neigs < typemax(Int))
       @warn "We found fewer eigvals than requested. Try increasing domain, or decreasing `tol`. This warning can be disabled with `sanity_check=false`." S
    end

    return (λgood,Vgood)
end

# Trapezoidal rule for a periodic function f
function ptrapz(f,g,a,b,N,pg)
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S0 = zero(f(t[1])); S1 = zero(S0)
    for i = 1:N
        temp = f(t[i])
        S0 += temp
        S1 += temp*g(t[i])
        next!(pg)
    end
    return h*S0, h*S1;
end


function ptrapz_parallel(f,g,a,b,N,channel)
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S = @distributed (+) for i = 1:N
        temp = f(t[i])
        put!(channel,true)
        [temp,temp*g(t[i])]
    end
    return h*S[1], h*S[2];
end
