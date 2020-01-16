using Test
using Iris

@testset "Common" begin
    include("Common/points.jl")
    include("Common/vectorfields.jl")
    include("Common/shapes.jl")
    include("Common/boundarylayers.jl")
    include("Common/boundaryconditions.jl")
    include("Common/boundaries.jl")
    include("Common/dielectricfunctions.jl")
    include("Common/pumpfunctions.jl")
    include("Common/lattices.jl")
    include("Common/dispersions.jl")
    include("Common/domains.jl")
    include("Common/laplacians.jl")
    # include("Common/curlcurls.jl")
    include("Common/selfenergies.jl")
    include("Common/simulations.jl")
    include("Common/helmholtzoperators.jl")
    include("Common/lufactorizations.jl")
end

@testset "Spectral" begin
end

@testset "Floquet" begin
end

@testset "Scattering" begin
end

@testset "SMatrices" begin
end

@testset "Lasing" begin
end

@testset "SaturableCPA" begin
end

@testset "TimeDomain" begin
end
