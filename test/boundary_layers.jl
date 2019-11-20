@testset "BoundaryLayers.jl" begin
    bl1 = noBL{1}(.1)
    bl2 = PML{2}(.1)
    bl3 = cPML{3}(.1)
    @test Iris.Common.BoundaryLayers.get_side(bl1)==1
    @test Iris.Common.BoundaryLayers.get_side(bl2)==2
    @test Iris.Common.BoundaryLayers.get_side(bl3)==3
    @test bl1(Point(.05))==0
    @test bl2(Point(.05))==0
end
