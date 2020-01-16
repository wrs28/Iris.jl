using Iris
using Test

@testset "BoundaryConditions" begin
    bc = noBC{1}()
    @test BCHermiticity(bc)<:HermitianBC
    @test Iris.Common.BCLocality(bc)<:Iris.Common.LocalBC
    @test Iris.Common.getside(bc) == 1
    bc = DirichletBC{2}()
    @test BCHermiticity(bc)<:HermitianBC
    @test Iris.Common.BCLocality(bc)<:Iris.Common.LocalBC
    @test Iris.Common.getside(bc) == 2
    bc = NeumannBC{3}()
    @test BCHermiticity(bc)<:HermitianBC
    @test Iris.Common.BCLocality(bc)<:Iris.Common.LocalBC
    @test Iris.Common.getside(bc) == 3
    bc = FloquetBC{4}()
    @test BCHermiticity(bc)<:HermitianBC
    @test Iris.Common.BCLocality(bc)<:Iris.Common.NonLocalBC
    @test Iris.Common.getside(bc) == 4
    bc = MatchedBC{5}(in=[1])
    @test BCHermiticity(bc)<:NonHermitianBC
    @test Iris.Common.BCLocality(bc)<:Iris.Common.NonLocalBC
    @test Iris.Common.getside(bc) == 5
    conj(bc).out == bc.in
end
