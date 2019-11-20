@testset "Shapes.jl" begin
    i = Interval(0,1)
    p = Point(.5)
    @test i(p)
    @test !i(3)
    @test i.origin==p
end
