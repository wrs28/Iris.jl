using Test
using Iris

@testset "Common" begin
    include("common.jl")
end

@testset "Spectral" begin
    include("spectral.jl")
end

@testset "Floquet" begin
    include("floquet.jl")
end

@testset "Scattering" begin
    include("scattering.jl")
end

@testset "SMatrices" begin
    include("smatrices.jl")
end

@testset "Lasing" begin
    include("lasing.jl")
end

@testset "SaturableCPA" begin
    include("saturablecpa.jl")
end

@testset "TimeDomain" begin
    include("timedomain.jl")
end
