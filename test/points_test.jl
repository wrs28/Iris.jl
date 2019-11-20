@testset "Points.jl" begin
    p = Point(1,2,3,4)
    @test 3*p == Point(3,6,9,12)
    @test p*3 == Point(3,6,9,12)
    @test p+p==2p
    @test p-p==0p
    @test ndims(p)==4
    @test p[2]==2
    @test firstindex(p)==1
    @test p[end]==4
    @test length(p)==4
    @test size(p)==(4,)
    @test eltype(p)==Float64
    @test p.x==1
    @test p.y==2
    @test p.z==3
    @test p.r==hypot(1,2,3,4)
    p = Point(3,4)
    @test p.r==5
    @test p.φ≈0.92729521
end
