# [Building a Simulation](@id common_manual)

A simulation comprises a list of domains. There are three kinds of domains: [`LatticeDomain`](@ref), [`NondispersiveDomain`](@ref), [`DispersiveDomain`](@ref).

A [`LatticeDomain`](@ref) defines a region that shares a common lattice. Therefore
a lattice domain contains the geometric information, such as position, about the sites.
It also defines the "background" index of refraction (e.g. for a scatterer in free space, n=1).
A simulation having only a single [`LatticeDomain`](Common.Domains.@ref) will be [`Symmetric`](@ref), in the sense
that the differential Maxwell operator will be symmetric.
At present, [`Unsymmetric`](@ref) simulations (i.e. those with two or more [`LatticeDomain`](@ref)s,
are only defined in dimension 2 or greater).

A [`NondispersiveDomain`](@ref) defines a region that shares a common dielectric function.
Typically this will be a piecewise-constant function, but can also be a user-defined function (see [Delectric Function](@ref) below).

A [`DispersiveDomain`](@ref) defines a region that shares a common (nonlinear) dispersion.
At the moment, only [two-level dispersions](@ref TwoLevelSystems) have been defined.

These regions do not need to be mutually exclusive. A given site will be assigned
the dielectric of its associated [`NondispersiveDomain`](@ref), defaulting to that defined
in its associated [`LatticeDomain`](@ref), with an additional dispersive (nonlinear) susceptability
defined by the associated [`DispersiveDomain`](@ref), defaulting to no additional
susceptability if the site is not contained within any [`DispersiveDomain`](@ref)s.
Ambiguities are resolved by assigning the domain that appears first in the arguments
provided to the [`Simulation`](@ref) constructor.

Once a series of domains has been defined, the constructor `Simulation(domains...) -> sim` fully
a simulation, which can then be passed to any of the main functions described in [`Spectral Analysis`](@ref spectral_manual), [`Band Dispersions`](@ref floquet_manual),
[`Scattering`](@ref scattering_manual), [`S-Matrix`](@ref smatrices_manual), [`Lasing`](@ref lasing_manual), [`Saturable CPA`](@ref saturablecpa_manual), [`FDTD`](@ref timedomain_manual).

## [Building a Domain](@id build_domain)
### [Lattice Domain](@ref Main.Iris.Common.Domains.LatticeDomain)
#### Boundary
#### Lattice
### [Nondispersive Domain](@ref NondispersiveDomain)
#### Shape
#### Dielectric Function
### [Dispersive Domain](@ref DispersiveDomain)
#### Shape
#### Dispersive Domain
