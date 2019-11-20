@testset "ElectricFields.jl" begin
    N = 10
    M = 3
    r = rand(ComplexF64,3N,M)
    pos = [Point(rand(Int,N)...) for i ∈ 1:N]
    e = ElectricField(pos,r)
    @test e==conj(conj(e))
    @test e + e == 2e
    @test e - e == 0e
    @test 4e/2==2e
    @test firstindex(e)==1
    @test lastindex(e)==3M*N
    @test length(e) == 3M*N
    @test size(e) == (3N,M)
    @test eltype(e) == ComplexF64
    @test typeof(e(2))<:ElectricField

end
