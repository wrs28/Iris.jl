using Iris
using Test

@testset "Self Energies" begin
    @testset "1D Symmetric" begin
        shape = Interval(-1,1)
        bnd = Boundary(shape,DirichletBC)
        lat = Lattices.Cartesian(.1)
        dom = LatticeDomain(bnd,lat)
        Σ = @test_nowarn Iris.Common.SelfEnergies.SelfEnergy{Iris.Symmetric}(dom,ones(ComplexF64,length(dom.x)+1))
        for i ∈ 1:20, j ∈ 1:20
            if i==1 && j==1
                @test Σ.Σ0[1][i,j] ≈ 100
                @test Σ.Σ1[1][1,1] ≈ 10
            elseif i==20 && j==20
                @test Σ.Σ0[2][20,20] ≈ 100
                @test Σ.Σ1[2][20,20] ≈ -10
            else
                @test Σ.Σ1[1][i,j] ≈ 0
                @test Σ.Σ1[2][i,j] ≈ 0
            end
        end
    end
    @testset "2D Symmetric" begin
        shape = Square(1,0,0)
        bnd = Boundary(shape,MatchedBC{1}(in=[1]),MatchedBC{2}(out=[1]),DirichletBC{3}(),DirichletBC{4}(),PML;depth=.3)
        lat = Lattices.Cartesian(.05,.05)
        dom = LatticeDomain(bnd,lat)
        s = @test_nowarn Iris.Common.SelfEnergies.SelfEnergy{Iris.Symmetric}(dom,ones(ComplexF64,20),ones(ComplexF64,20))
    end
end
