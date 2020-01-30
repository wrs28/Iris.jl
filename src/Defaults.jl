# file contains default parameters used throughout Iris

# PML/cPML parameters
const EXTINCTION = 1e-10 # extinction in PML layer
const SCALING_ANGLE = .15 # phase in conductivity to accelerate evanscent decay
const BL_DEPTH = .1 # default depth of cPML or PML layer
const DEFAULT_BC = :noBC # default boundary condition (be specified as symbol or string)
const DEFAULT_BL = :noBL # default boundary layer (specified as symbol or string)
const LUPACK = :USolver # default LU package, must be one of USolver, MSolver, PSolver (specified as symbol or string)
const INDEX_OFFSET = 7 # position index offset from midpoint where field is real (SPECTRAL), or real and non-vanishing (LASING/SCPA)
const EQUIVALENT_SOURCE_RELATIVE_CUTOFF = 1e-8 # relative tolerance condition for dropping components of equivalent current
const DEFAULT_CFL_NUMBER = 0.9
const DEFAULT_LINEAR_EIGENSOLVER = :Arpack
const DEFAULT_NONLINEAR_EIGENSOLVER = :NonlinearEigenproblems
const NUM_SUBPIXELS = 5
const ORTHOGONALIZE_OVERLAP_THRESHOLD = .1

# 2D unsymmetric shape defaults
const NL_NORMAL_ALGORITHM = :LD_MMA # NLopt algorithm used to find the normal, tangent, and distance to shape perimeter
const NL_NORMAL_XTOL = 1e-5 # NLopt's xtol_rel
const NL_NORMAL_FTOL = 1e-5 # NLopt's ftol_rel
const NL_NORMAL_MAXEVAL = 400 # NLopt's maxeval

# for plotting
const SHAPE_FILL_ALPHA = .2
const SHAPE_COLOR = :red
const BOUNDARY_COLOR = :gray
const BOUNDARY_PML_COLOR = :blue
const BOUNDARY_CPML_COLOR = :red
const BOUNDARY_DIRICHLET_LINETYPE = :solid
const BOUNDARY_NEUMANN_LINETYPE = :dash
const BOUNDARY_FLOQUET_MATCHED_LINETYPE = :dashdot
const BOUNDARY_BC_COLOR = :black
const MARKERSHAPE = :rect
const LINK_SCALE_REDUCTION = 1/2
const MARKERSIZE_SCALE = 135
const PLOT_QUANTILE = .95 # optimize the dynamic range of plotting wavefunctions to cover this percentile
const PLOT_SCALE_FUDGE = 1.3 # multiply plot_quantile result by this factor to get colorbar scale

# Pretty Printing Colors
const PRINTED_COLOR_NUMBER = :light_cyan # numbers
const PRINTED_COLOR_VARIABLE = :cyan # varialbes and functions
const PRINTED_COLOR_WARN = :light_yellow # to grab attention but not necessarily b/c of an error
const PRINTED_COLOR_GOOD = :green # successful solution
const PRINTED_COLOR_DARK = 63 # names of structures that are part of simulation, but do not themselves contain calculations
const PRINTED_COLOR_LIGHT = 171 # names of structures used to calculate and/or containing solutions
const PRINTED_COLOR_BAD = :light_red # indicates error
const PRINTED_COLOR_INSTRUCTION = 244 # parenthetical instructions

"""
PML extinction factor: `$EXTINCTION`
"""
EXTINCTION

"""
PML scaling angle: `$SCALING_ANGLE`
"""
SCALING_ANGLE

"""
default PML depth: `$BL_DEPTH`
"""
BL_DEPTH

"""
default boundary condition: `$DEFAULT_BC`
"""
DEFAULT_BC

"""
default boundary layer: `$DEFAULT_BL`
"""
DEFAULT_BL

"""
default LU package (one of `USolver`, `PSolver`, `MSolver` for UMFPACK, Pardiso, MUMPS): `$LUPACK`
"""
LUPACK

"""
index offset from midpoint where field is real (SPECTRAL), or real and nonzero (LASING/SaturableCPA): `$INDEX_OFFSET`
"""
INDEX_OFFSET

"""
relative tolerance condition for dropping components of equivalent current: `$EQUIVALENT_SOURCE_RELATIVE_CUTOFF`
"""
EQUIVALENT_SOURCE_RELATIVE_CUTOFF

"""
Courant-Friedrichs-Lewy number (N*space-to-time step ratio in time domain simulations in `N` dimensions): `$DEFAULT_CFL_NUMBER`
"""
DEFAULT_CFL_NUMBER

"""
`$DEFAULT_LINEAR_EIGENSOLVER`
"""
DEFAULT_LINEAR_EIGENSOLVER

"""
`$DEFAULT_NONLINEAR_EIGENSOLVER`
"""
DEFAULT_NONLINEAR_EIGENSOLVER

"""
subsampling rate for simulation smoothing: `$NUM_SUBPIXELS`
"""
NUM_SUBPIXELS

#
# # Tessellation defaults
# const TESSELLATION_EXPANSION_FACTOR = 1.025 # to make accessing the tessellation a little more stable
#
# # for interpolating in Simulations
# const N_CUBIC = 10
# const N_QUADRATIC = 6
# const N_BILINEAR = 4
# const N_LINEAR = 3
#
# const NNN_BULK = N_CUBIC # number of nearest neighbors used in bulk
# const NNN_SURFACE = N_LINEAR # number of nearest neighbors for surface
# const NNN_CORNER = N_LINEAR # number of nearest neighbors at corners
