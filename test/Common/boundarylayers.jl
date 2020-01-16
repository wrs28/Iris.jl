using Iris
using Test

@testset "BoundaryLayers" begin
    pml1 = PML{1}(.1)
    @test_nowarn pml1.start = 0
    @test_nowarn pml1.stop = pml1.start - pml1.depth
    pml2 = PML{2}(.1)
    @test_nowarn pml2.start = 0
    @test_nowarn pml2.stop = pml2.start + pml2.depth
    @test pml1(-.5) ≈ pml2(.5)
    @test Iris.Common.getside(pml1) == 1
    @test Iris.Common.getside(pml2) == 2

    cpml1 = cPML{1}(.1)
    @test_nowarn cpml1.start = 0
    @test_nowarn cpml1.stop = cpml1.start - cpml1.depth
    cpml2 = cPML{2}(.1)
    @test_nowarn cpml2.start = 0
    @test_nowarn cpml2.stop = cpml2.start + cpml2.depth
    @test cpml1(-.5) ≈ cpml2(.5)
    @test Iris.Common.getside(cpml1) == 1
    @test Iris.Common.getside(cpml2) == 2

    @test pml1(-.5) ≈ -conj(cpml2(.5))

    nbl1 = noBL{1}()
    @test_nowarn nbl1.start = 0
    @test_nowarn nbl1.stop = nbl1.start - nbl1.depth
    nbl2 = noBL{2}(.1)
    @test_nowarn nbl2.start = 0
    @test_nowarn nbl2.stop = nbl2.start + nbl2.depth
    @test nbl1(-.5) ≈ conj(nbl2(-.5))
    @test Iris.Common.getside(nbl1) == 1
    @test Iris.Common.getside(nbl2) == 2
end
