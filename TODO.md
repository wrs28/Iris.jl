# To Do

* [x] fix αε in PML
* [x] Get working eig_kl
* [x] Get working eig_knl
* [x] Get working eig_cf
* [x] check PML normal incidence
* [x] check PML nonnormal incidence
* [x] boundary matching
* [x] relearn NLsolve, see if still best
* [ ] disperions: reorganize, two-level, kerr/chi-3
* [x] build single-mode SALT
* [x] define SCPA through SALT
* [ ] floquet/periodic

In no particular order

* [ ] Add periodic (floquet) boundary conditions
* [ ] Band structure calculations
* [x] SALT
* [x] Saturable CPA
* [ ] Polar Coordinates
* [ ] Nearest Neighbor shortcut in polar coordinates
* [x] Boundary Layers
* [x] Boundary matching
* [x] Make nonlinear solvers use Pardiso/MUMPS
* [ ] Fix nonlinear contour solver with branch cuts, seek out only physical solutions

For tonight
* [x] Nonlinear Scattering: nonlinear scattering struct, use nlsolve from there
* [x] Revisit multiple RHS/LHS for ArnoldiMethodTransformations in UMFPACK, Pardiso, MUMPS cases
* [x] Reduce number of allocations in ArnoldiMethodTransformations
* [x] Reduce number of allocations in maxwell_lep, maxwell
* [ ] Jacobian for nonlinear scattering
* [ ] Jacobian for SALT/SCPA
* [x] Plots for: Electric Field, Simulation, Scattering solution
* [ ] Threshold finder for SALT
* [ ] Bootstrap for SALT
* [ ] Time Domain for Maxwell-Bloch (5-component field)
