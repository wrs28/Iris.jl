using Iris
using Test

@testset "Boundaries" begin
    @testset "1D" begin
        shape = Interval(0,1)
        @test_nowarn Boundary(shape,DirichletBC)
        @test_nowarn Boundary(shape,PML)
        @test_nowarn Boundary(shape,DirichletBC,PML{1}(.1),noBL{2}())
        bcs = (DirichletBC{1}(),DirichletBC{2}())
        bls = (PML{1}(.1),PML{1}(.1),PML{1}(.1))
        @test_throws ArgumentError Boundary(shape,bcs,bls)
        bls = (noBL{1}(.1),PML{2}(.2))
        bcs = (DirichletBC{1}(),DirichletBC{2}(),DirichletBC{3}())
        @test_throws ArgumentError Boundary(shape,bcs,bls)
        bcs = (DirichletBC{1}(),DirichletBC{2}())
        @test_nowarn bnd = Boundary(shape,bcs,bls)
        bnd = Boundary(shape,bcs,bls)
        @test ndims(bnd)==1
        @test Iris.Common.Boundaries.nsides(bnd)==2
    end
    @testset "2D" begin
        @testset "Square" begin
            shape = Square(1,0,0)
            @test_nowarn Boundary(shape,DirichletBC)
            @test_nowarn Boundary(shape,PML)
            @test_nowarn Boundary(shape,DirichletBC,PML{1}(.1),noBL{2}())
            bcs = (DirichletBC{1}(),DirichletBC{2}())
            bls = (PML{1}(),PML{1}(.1),PML{1}(.1),cPML{2}(.3),noBL{4}())
            @test_throws ArgumentError Boundary(shape,bcs,bls)
            bls = (noBL{1}(.1),PML{2}(.2))
            bcs = (DirichletBC{1}(),DirichletBC{2}(),DirichletBC{3}(),MatchedBC{3}(),NeumannBC{1}())
            @test_throws ArgumentError Boundary(shape,bcs,bls)
            bls = (noBL{3}(.1),PML{4}(.2),cPML{2}(),noBL{1}())
            bcs = (DirichletBC{2}(),NeumannBC{1}(),DirichletBC{4}(),DirichletBC{3}())
            @test_nowarn bnd = Boundary(shape,bcs,bls)
            bnd = Boundary(shape,bcs,bls)
            @test ndims(bnd)==2
            @test Iris.Common.Boundaries.nsides(bnd)==4
        end
        @testset "Circle" begin
            shape = Circle(1,0,0)
            @test_nowarn Boundary(shape)
            @test_nowarn Boundary(shape,DirichletBC)
            @test_nowarn Boundary(shape,PML)
            @test_throws ArgumentError Boundary(shape,DirichletBC,PML{1}(.1),noBL{2}())
            bcs = (DirichletBC{1}(),)
            bls = (PML{1}(),)
            bnd = Boundary(shape,bcs,bls)
            @test_nowarn Boundary(shape,bcs,bls)
            bls = (noBL{1}(.1),)
            bcs = (DirichletBC{1}(),)
            @test_nowarn Boundary(shape,bcs,bls)
            bls = (noBL{3}(.1),PML{4}(.2),cPML{2}(),noBL{1}())
            bcs = (DirichletBC{2}(),NeumannBC{1}(),DirichletBC{4}(),DirichletBC{3}())
            @test_throws ArgumentError bnd = Boundary(shape,bcs,bls)
            @test ndims(bnd)==2
            @test Iris.Common.Boundaries.nsides(bnd)==1
        end
    end
end
