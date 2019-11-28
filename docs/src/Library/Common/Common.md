# Common

Core routines for creating a `Simulation` object.
Many default parameters defined as [Global Defaults](@ref).

See [Building a Simulation](@ref) for tutorial.

---

Exported methods:
  * [`Point`](@ref)
  * [`ElectricField`](@ref)
  * [`AbstractShape`](@ref)
  * [`Interval`](@ref)
  * [`PML`](@ref)
  * [`cPML`](@ref)
  * [`noBL`](@ref)
  * [`DirichletBC`](@ref)
  * [`NeumannBC`](@ref)
  * [`FloquetBC`](@ref)
  * [`MatchedBC`](@ref)
  * [`noBC`](@ref)
  * [`BCHermiticity`](@ref)
  * [`HermitianBC`](@ref)
  * [`NonHermitianBC`](@ref)
  * [`Boundary`](@ref)
  * [`AbstractDielectricFunction`](@ref)
  * [`AbstractPumpFunction`](@ref)
  * [`TwoLevelSystem`](@ref)
  * [`susceptability`](@ref)
  * [`Lattices`](@ref)
  * [`Lattice`](@ref)
  * [`LatticeDomain`](@ref)
  * [`NondispersiveDomain`](@ref)
  * [`DispersiveDomain`](@ref)
  * [`Curlcurl`](@ref)
  * [`Simulation`](@ref)
  * [`smooth!`](@ref)
  * [`update!`](@ref)
  * [`update_dielectric!`](@ref)
  * [`update_pump!`](@ref)
  * [`Unsymmetric`](@ref)
  * [`Symmetric`](@ref)
  * [`Hermitian`](@ref)
  * [`Maxwell`](@ref)
  * [`DEFAULT_LUPACK`](@ref)
  * [`MSolver`](@ref)
  * [`PSolver`](@ref)
  * [`USolver`](@ref)
  * [`AbstractLUPACK`](@ref)

---

```@docs
Common
```
---

### Global Defaults
```@autodocs
Modules = [Common]
Order = [:constant]
```
