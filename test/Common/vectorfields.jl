using Iris
using Test

@testset "VectorFields" begin
    @testset "VectorField" begin
        N = 10 # dimension
        M = 2 # number of components
        n = 15 # number of sites
        nmodes = 3
        r = rand(ComplexF64,M*n,nmodes)
        pos = Vector{Point{N}}(undef,n)
        for i ∈ 1:n pos[i] = Point(randn(Float64,N)...) end
        start = Point(zeros(Int,N)...)
        stop = Point(ones(Int,N)...)
        start_inds = rand(Int,N)
        stop_inds = rand(Int,N)
        e = VectorField{M}(pos,r,start,stop,start_inds,stop_inds)
        @test e==conj(conj(e))
        @test e + e == 2e
        @test e - e == 0e
        @test 4e/2==2e
        @test firstindex(e)==1
        @test lastindex(e)==M*n*nmodes
        @test length(e) == M*n*nmodes
        @test size(e) == (M*n,nmodes)
        @test eltype(e) == ComplexF64
        @test typeof(e(2))<:VectorField
        @test_nowarn VectorField(e,rand(ComplexF64,size(e)...))
        @test_nowarn update!(e,rand(ComplexF64,size(e)...))
        @test_nowarn component(1,e)
        @test_throws BoundsError component(3,e; verbose=false)
    end
    @testset "ScalarField" begin
        N = 1 # dimension
        M = 1 # number of components
        n = 25 # number of sites
        nmodes = 4
        r = rand(ComplexF64,M*n,nmodes)
        pos = Vector{Point{N}}(undef,n)
        for i ∈ 1:n pos[i] = Point(randn(Float64,N)...) end
        start = Point(zeros(Int,N)...)
        stop = Point(ones(Int,N)...)
        start_inds = rand(Int,N)
        stop_inds = rand(Int,N)
        e = ScalarField(pos,r,start,stop,start_inds,stop_inds)
        @test e==conj(conj(e))
        @test e + e == 2e
        @test e - e == 0e
        @test 4e/2==2e
        @test firstindex(e)==1
        @test lastindex(e)==M*n*nmodes
        @test length(e) == M*n*nmodes
        @test size(e) == (M*n,nmodes)
        @test eltype(e) == ComplexF64
        @test typeof(e(2))<:ScalarField
        @test_nowarn ScalarField(e,rand(ComplexF64,size(e)...))
        @test_nowarn update!(e,rand(ComplexF64,size(e)...))
        @test_nowarn component(1,e)
        @test_throws BoundsError component(3,e; verbose=false)
    end
    @testset "ElectricField" begin
        N = 1 # dimension
        M = 3 # number of components
        n = 25 # number of sites
        nmodes = 4
        r = rand(ComplexF64,M*n,nmodes)
        pos = Vector{Point{N}}(undef,n)
        for i ∈ 1:n pos[i] = Point(randn(Float64,N)...) end
        start = Point(zeros(Int,N)...)
        stop = Point(ones(Int,N)...)
        start_inds = rand(Int,N)
        stop_inds = rand(Int,N)
        e = ElectricField(pos,r,start,stop,start_inds,stop_inds)
        @test e==conj(conj(e))
        @test e + e == 2e
        @test e - e == 0e
        @test 4e/2==2e
        @test firstindex(e)==1
        @test lastindex(e)==M*n*nmodes
        @test length(e) == M*n*nmodes
        @test size(e) == (M*n,nmodes)
        @test eltype(e) == ComplexF64
        @test typeof(e(2))<:ElectricField
        @test_nowarn ElectricField(e,rand(ComplexF64,size(e)...))
        @test_nowarn update!(e,rand(ComplexF64,size(e)...))
        @test_nowarn component(1,e)
        @test_throws BoundsError component(4,e; verbose=false)
    end
end
