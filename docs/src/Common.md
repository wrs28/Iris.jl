# Building a Simulation

A simulation comprises a list of domains. There are three kinds of domains: a [`LatticeDomain`](@ref), [`NondispersiveDomain`](@ref), [`DispersiveDomain`](@ref).

A [`LatticeDomain`](@ref) defines a region that shares a common lattice.

## 1-Dimensional

The [`LatticeDomain`](@ref) in 1-D has the following signature:

````JULIA
LatticeDomain(::Boundary{1}, ::Lattice{1}, [n=1; type, name, fit=true) -> domain
````
