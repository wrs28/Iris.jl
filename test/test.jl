using Iris
using SpecialFunctions

c = Circle(1,2,3)
lat = Lattices.Cartesian(1/21,1/21)
bnd = Boundary(c,DirichletBC())
sim = Simulation(Cavity(bnd,lat))
k,ψ = eig_kl(sim,eps(),nev=20)
@testset "Zeros of Bessel J" begin
    @test abs(besselj(0,k[1])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(1,k[2])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(1,k[3])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(2,k[4])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(2,k[5])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(0,k[6])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(3,k[7])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(3,k[8])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(1,k[9])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(1,k[10])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(4,k[11])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(4,k[12])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(2,k[13])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(2,k[14])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(0,k[15])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(5,k[16])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(5,k[17])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(3,k[18])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(3,k[19])) ≤ max(lat.dx,lat.dy)
    @test abs(besselj(6,k[20])) ≤ max(lat.dx,lat.dy)
end

r = Rectangle(1,2,3,4,.1)
bnd = Boundary(r,NeumannBC())
sim = Simulation(Cavity(bnd,lat))
k,ψ = eig_kl(sim,eps(),nev=5)
@testset "1:2 Rectangle" begin
    @test abs((k[1])^2 - (0π/1)^2 - (0π/2)^2) ≤ 100*2*max(lat.dx,lat.dy)^2
    @test abs((k[2])^2 - (0π/1)^2 - (1π/2)^2) ≤ 100*2*max(lat.dx,lat.dy)^2
    @test abs((k[3])^2 - (0π/1)^2 - (2π/2)^2) ≤ 100*2*max(lat.dx,lat.dy)^2
    @test abs((k[4])^2 - (1π/1)^2 - (0π/2)^2) ≤ 100*2*max(lat.dx,lat.dy)^2
    @test abs((k[5])^2 - (1π/1)^2 - (1π/2)^2) ≤ 100*2*max(lat.dx,lat.dy)^2
end
