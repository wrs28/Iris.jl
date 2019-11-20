"""
    module TimeDomain
"""
module TimeDomain

export maxwell_bloch

include("Curls.jl")
using .Curls

include("PolarizationFields.jl")
using .PolarizationFields

include("InversionFields.jl")
using .InversionFields

# include("AuxilliaryFields.jl")
# using .AuxilliaryFields

# include("Maxwell/Maxwell.jl")
# using .MaxwellFDTD

include("MaxwellBloch/MaxwellBloch.jl")
using .MaxwellBlochFDTD

end # module
