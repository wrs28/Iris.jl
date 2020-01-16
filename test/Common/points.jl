using Iris
using Test

@testset "Points" begin
    @testset "1D Cartesian" begin
        p = Point(1)
        @test_throws ArgumentError Point{Cartesian}(Inf)
        @test_throws ArgumentError Point{Cartesian}(NaN)
        @test 3*p == Point{Cartesian}(3)
        @test p*3 == Point{Cartesian}(3)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==1
        @test p[1] ≈ 1
        @test firstindex(p)==1
        @test p[end] ≈ 1
        @test length(p)==1
        @test size(p)==(1,)
        @test eltype(p)==Float64
        @test p.x ≈ 1
        @test p.y ≈ 0
        @test p.z ≈ 0
        @test p.s ≈ 1
        @test p.r ≈ hypot(1)
        @test p.ϕ ≈ 0
        @test p.θ ≈ π/2
    end
    @testset "2D Cartesian" begin
        p = Point(1,2.)
        @test_throws ArgumentError Point{Cartesian}(Inf,2.)
        @test_throws ArgumentError Point{Cartesian}(1,NaN)
        @test 3*p == Point{Cartesian}(3,6)
        @test p*3 == Point{Cartesian}(3,6)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==2
        @test p[2]==2
        @test firstindex(p)==1
        @test p[end]==2
        @test length(p)==2
        @test size(p)==(2,)
        @test eltype(p)==Float64
        @test p.x==1
        @test p.y==2
        @test p.r==hypot(1,2)

        p = Point(3,4)
        @test p.x ≈ 3
        @test p.y ≈ 4
        @test p.z ≈ 0
        @test p.s ≈ 5
        @test p.r ≈ 5
        @test p.ϕ ≈ 0.92729521
        @test p.θ ≈ π/2
    end
    @testset "3D Cartesian" begin
        p = Point(1,2,3.)
        @test_throws ArgumentError Point{Cartesian}(Inf,2.,3)
        @test_throws ArgumentError Point{Cartesian}(1,NaN,3)
        @test 3*p == Point{Cartesian}(3,6,9)
        @test p*3 == Point{Cartesian}(3,6,9)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==3
        @test p[2]==2
        @test firstindex(p)==1
        @test p[end]==3
        @test length(p)==3
        @test size(p)==(3,)
        @test eltype(p)==Float64

        p = Point{Cartesian}(1,sqrt(8),4)
        @test p.x ≈ 1
        @test p.y ≈ sqrt(8)
        @test p.z ≈ 4
        @test p.s ≈ 3
        @test p.r ≈ 5
        @test p.ϕ ≈ 1.2309594173407747
        @test p.θ ≈ 0.6435011087932843
    end
    @testset "4D Cartesian" begin
        p = Point(1,2,3,4.)
        @test_throws ArgumentError Point{Cartesian}(Inf,2.,3,4.)
        @test_throws ArgumentError Point{Cartesian}(1,NaN,3,4.)
        @test 3*p == Point{Cartesian}(3,6,9,12)
        @test p*3 == Point{Cartesian}(3,6,9,12)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==4
        @test p[2]==2
        @test firstindex(p)==1
        @test p[end]==4
        @test length(p)==4
        @test size(p)==(4,)
        @test eltype(p)==Float64

        p = Point{Cartesian}(1,sqrt(8),4,12)
        @test p.x ≈ 1
        @test p.y ≈ sqrt(8)
        @test p.z ≈ 4
        @test p.s ≈ 3
        @test p.r ≈ 13
        @test_throws ErrorException p.ϕ
        @test_throws ErrorException p.θ
    end
    @testset "1D Polar" begin
        p = Point{Polar}(1)
        @test_throws ArgumentError Point{Polar}(Inf)
        @test_throws ArgumentError Point{Polar}(NaN)
        @test_throws ArgumentError Polar(-1)
        @test 3*p == Point{Polar}(3)
        @test p*3 == Point{Polar}(3)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==1
        @test firstindex(p)==1
        @test p[end]==1
        @test length(p)==1
        @test size(p)==(1,)
        @test eltype(p)==Float64

        p = Point{Polar}(1)
        @test p.x ≈ 1
        @test p.y ≈ 0
        @test p.z ≈ 0
        @test p.s ≈ 1
        @test p.r ≈ 1
        @test p.ϕ ≈ 0
        @test p.θ ≈ π/2
    end
    @testset "2D Polar" begin
        p = Point{Polar}(1,2.)
        @test_throws ArgumentError Point{Polar}(Inf,2.)
        @test_throws ArgumentError Point{Polar}(NaN,1)
        @test_throws ArgumentError Polar(-1,2.)
        @test 3*p == Point{Polar}(3,2)
        @test p*3 == Point{Polar}(3,2)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==2
        @test p[2]==2
        @test firstindex(p)==1
        @test p[end]==2
        @test length(p)==2
        @test size(p)==(2,)
        @test eltype(p)==Float64

        p = Point{Polar}(1,π/3)
        @test p.x ≈ 1/2
        @test p.y ≈ sqrt(3)/2
        @test p.z ≈ 0
        @test p.s ≈ 1
        @test p.r ≈ 1
        @test p.ϕ ≈ π/3
        @test p.θ ≈ π/2
    end
    @testset "3D Polar" begin
        p = Point{Polar}(1,2,3.)
        @test_throws ArgumentError Point{Polar}(Inf,2.,NaN)
        @test_throws ArgumentError Point{Polar}(NaN,1,3.)
        @test_throws ArgumentError Point{Polar}(-1,2.,3.)
        @test 3*p == Point{Polar}(3,2,9)
        @test p*3 == Point{Polar}(3,2,9)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==3
        @test p[2]==2
        @test firstindex(p)==1
        @test p[end]==3
        @test length(p)==3
        @test size(p)==(3,)
        @test eltype(p)==Float64

        p = Point{Polar}(3,π/3,4)
        @test p.x ≈ 3/2
        @test p.y ≈ 3*sqrt(3)/2
        @test p.z ≈ 4
        @test p.s ≈ 3
        @test p.r ≈ 5
        @test p.ϕ ≈ π/3
        @test p.θ ≈ 0.6435011087932843
    end
    @testset "1D Spherical" begin
        p = Point{Spherical}(1)
        @test_throws ArgumentError Point{Spherical}(NaN)
        @test_throws ArgumentError Point{Spherical}(Inf)
        @test_throws ArgumentError Point{Spherical}(-1)
        @test 3*p ≈ Point{Spherical}(3)
        @test p*3 ≈ Point{Spherical}(3)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==1
        @test firstindex(p)==1
        @test p[end]==1
        @test length(p)==1
        @test size(p)==(1,)
        @test eltype(p)==Float64

        p = Point{Spherical}(3)
        @test p.x ≈ 3
        @test p.y ≈ 0
        @test p.z ≈ 0
        @test p.s ≈ 3
        @test p.r ≈ 3
        @test p.ϕ ≈ 0
        @test p.θ ≈ 0
    end
    @testset "2D Spherical" begin
        p = Point{Spherical}(1,π/3)
        @test_throws ArgumentError Point{Spherical}(NaN,1)
        @test_throws ArgumentError Point{Spherical}(Inf,1)
        @test_throws ArgumentError Point{Spherical}(-1,1)
        @test 3*p ≈ Point{Spherical}(3,π/3)
        @test p*3 ≈ Point{Spherical}(3,π/3)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==2
        @test firstindex(p)==1
        @test p[end]==π/3
        @test length(p)==2
        @test size(p)==(2,)
        @test eltype(p)==Float64

        p = Point{Spherical}(3,π/3)
        @test p.x ≈ 3*sqrt(3)/2
        @test p.y ≈ 0
        @test p.z ≈ 3/2
        @test p.s ≈ p.x
        @test p.r ≈ 3
        @test p.ϕ ≈ 0
        @test p.θ ≈ π/3
    end
    @testset "3D Spherical" begin
        p = Point{Spherical}(1,2,3)
        @test_throws ArgumentError Point{Spherical}(NaN,1,2)
        @test_throws ArgumentError Point{Spherical}(Inf,1,NaN)
        @test_throws ArgumentError Point{Spherical}(-1,1,1)
        @test 3*p == Point{Spherical}(3,2,3)
        @test p*3 == Point{Spherical}(3,2,3)
        @test p+p≈2p
        @test p-p≈0p
        @test ndims(p)==3
        @test p[2]==2
        @test firstindex(p)==1
        @test p[end]==3
        @test length(p)==3
        @test size(p)==(3,)
        @test eltype(p)==Float64

        p = Point{Spherical}(3,π/3,π/4)
        @test p.x ≈ 3*(sqrt(3)/2)*(1/sqrt(2))
        @test p.y ≈ 3*(sqrt(3)/2)*(1/sqrt(2))
        @test p.z ≈ 3/2
        @test p.s ≈ 3*sqrt(3)/2
        @test p.r ≈ 3
        @test p.ϕ ≈ π/4
        @test p.θ ≈ π/3
    end
end
