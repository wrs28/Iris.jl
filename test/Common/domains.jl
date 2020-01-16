using Iris
using Test

@testset "Domains" begin
    @testset "1D" begin
        shape = Interval(-1,1)
        bnd = Boundary(shape)
        @testset "Nondispersive" begin
            ε = @test_nowarn NondispersiveDomain(shape,1+.1im,:test,:only)
            @test_nowarn NondispersiveDomain(shape,1)
            @test_nowarn NondispersiveDomain(shape,1im)
            @test_nowarn NondispersiveDomain(shape,1,.1)
            @test_nowarn NondispersiveDomain(shape,1,0)
            @test_nowarn NondispersiveDomain(shape,0,1)
            @test_nowarn NondispersiveDomain(bnd)
        end
        @testset "Dispersive" begin
            disp = TwoLevelSystem(3,.1,1)
            @test_nowarn DispersiveDomain(bnd,disp;type=:test,name=:only)
            @test_nowarn DispersiveDomain(bnd,disp,.5;type=:test,name=:only)
            ε = @test_nowarn DispersiveDomain(shape,disp;type=:test,name=:only)
        end
        @testset "Lattice" begin
            @testset "Cartesian" begin
                lat = Lattices.Cartesian(.1)
                dom = @test_nowarn LatticeDomain(bnd,lat)
                @test_nowarn LatticeDomain(lat,bnd)
                @test_nowarn LatticeDomain(bnd,lat,2)
                @test_nowarn LatticeDomain(lat,bnd,2)
                @test_nowarn LatticeDomain(bnd,2,lat)
                @test_nowarn LatticeDomain(lat,2,bnd)
                @test_nowarn LatticeDomain(2,bnd,lat)
                @test_nowarn LatticeDomain(2,lat,bnd)
            end
            # @testset "Polar" begin
            #     lat = Lattices.Cartesian(.1)
            #     dom = @test_nowarn LatticeDomain(bnd,lat)
            #     @test_nowarn LatticeDomain(lat,bnd)
            #     @test_nowarn LatticeDomain(bnd,lat,2)
            #     @test_nowarn LatticeDomain(lat,bnd,2)
            #     @test_nowarn LatticeDomain(bnd,2,lat)
            #     @test_nowarn LatticeDomain(lat,2,bnd)
            #     @test_nowarn LatticeDomain(2,bnd,lat)
            #     @test_nowarn LatticeDomain(2,lat,bnd)
            # end
            # @testset "Spherical" begin
            #     lat = Lattices.Cartesian(.1)
            #     dom = @test_nowarn LatticeDomain(bnd,lat)
            #     @test_nowarn LatticeDomain(lat,bnd)
            #     @test_nowarn LatticeDomain(bnd,lat,2)
            #     @test_nowarn LatticeDomain(lat,bnd,2)
            #     @test_nowarn LatticeDomain(bnd,2,lat)
            #     @test_nowarn LatticeDomain(lat,2,bnd)
            #     @test_nowarn LatticeDomain(2,bnd,lat)
            #     @test_nowarn LatticeDomain(2,lat,bnd)
            # end
        end
    end
    @testset "2D" begin
        shape = Square(1,0,0)
        bnd = Boundary(shape)
        @testset "Nondispersive" begin
            ε = @test_nowarn NondispersiveDomain(shape,1+.1im,:test,:only)
            @test_nowarn NondispersiveDomain(shape,1)
            @test_nowarn NondispersiveDomain(shape,1im)
            @test_nowarn NondispersiveDomain(shape,1,.1)
            @test_nowarn NondispersiveDomain(shape,1,0)
            @test_nowarn NondispersiveDomain(shape,0,1)
            @test_nowarn NondispersiveDomain(bnd)
        end
        @testset "Dispersive" begin
            disp = TwoLevelSystem(3,.1,1)
            @test_nowarn DispersiveDomain(bnd,disp;type=:test,name=:only)
            @test_nowarn DispersiveDomain(bnd,disp,.5;type=:test,name=:only)
            ε = @test_nowarn DispersiveDomain(shape,disp;type=:test,name=:only)
        end
        @testset "Lattice" begin
            @testset "Cartesian" begin
                lat = Lattices.Cartesian(.1,.1)
                dom = @test_nowarn LatticeDomain(bnd,lat)
                @test_nowarn LatticeDomain(lat,bnd)
                @test_nowarn LatticeDomain(bnd,lat,2)
                @test_nowarn LatticeDomain(lat,bnd,2)
                @test_nowarn LatticeDomain(bnd,2,lat)
                @test_nowarn LatticeDomain(lat,2,bnd)
                @test_nowarn LatticeDomain(2,bnd,lat)
                @test_nowarn LatticeDomain(2,lat,bnd)
                lat = Lattices.Cartesian(.05,.05,angle=.05)
                @test_nowarn typeof(LatticeDomain(bnd,lat))<:LatticeDomain{2,Domains.Unsymmetric}
            end
            # @testset "Polar" begin
            #     lat = Lattices.Cartesian(.1)
            #     dom = @test_nowarn LatticeDomain(bnd,lat)
            #     @test_nowarn LatticeDomain(lat,bnd)
            #     @test_nowarn LatticeDomain(bnd,lat,2)
            #     @test_nowarn LatticeDomain(lat,bnd,2)
            #     @test_nowarn LatticeDomain(bnd,2,lat)
            #     @test_nowarn LatticeDomain(lat,2,bnd)
            #     @test_nowarn LatticeDomain(2,bnd,lat)
            #     @test_nowarn LatticeDomain(2,lat,bnd)
            # end
            # @testset "Spherical" begin
            #     lat = Lattices.Cartesian(.1)
            #     dom = @test_nowarn LatticeDomain(bnd,lat)
            #     @test_nowarn LatticeDomain(lat,bnd)
            #     @test_nowarn LatticeDomain(bnd,lat,2)
            #     @test_nowarn LatticeDomain(lat,bnd,2)
            #     @test_nowarn LatticeDomain(bnd,2,lat)
            #     @test_nowarn LatticeDomain(lat,2,bnd)
            #     @test_nowarn LatticeDomain(2,bnd,lat)
            #     @test_nowarn LatticeDomain(2,lat,bnd)
            # end
        end
    end
end
