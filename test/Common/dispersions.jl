using Iris
using Test

@testset "Dispersions" begin
    @testset "No Dispersion" begin
        @test_nowarn NoDispersion()
        d = NoDispersion()
        @test susceptability(d,1,2,3) ≈ 0
    end
    @testset "Two Level" begin
        @test_nowarn TwoLevelSystem()
        d = TwoLevelSystem()
        @test_nowarn susceptability(d,3)
        @test susceptability(d,3) ≈ 0
        ω = 3
        ωs = 2.5 .+ 1rand(3)
        @testset "Scalar" begin
            r = rand(ComplexF64,20,length(ωs))
            pos = Vector{Point{1}}(undef,20)
            for i ∈ 1:20 pos[i] = Point(randn(Float64,1)...) end
            e = ScalarField(pos,r)
            @test_nowarn susceptability(d,ω,ωs,e)
            d.D₀ = .1
            d.ωa = ω
            @test_nowarn susceptability(d,ω,ωs,e)
        end
        @testset "Vector" begin
            N = 1 # dimension
            M = 3 # number of components
            n = 25 # number of sites
            nmodes = 3
            r = rand(ComplexF64,M*n,nmodes)
            pos = Vector{Point{N}}(undef,n)
            for i ∈ 1:n pos[i] = Point(randn(Float64,N)...) end
            e = ElectricField(pos,r)
            @test_nowarn susceptability(d,ω,ωs,e)
            d.D₀ = .1
            d.ωa = ω
            @test_nowarn susceptability(d,ω,ωs,e)
        end
    end
    # @testset "Kerr" begin
    #     d = Kerr()
    #     @test susceptability(d,1,2,3) ≈ 0.0
    # end
end
