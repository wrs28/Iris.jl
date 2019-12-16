"""
Plot reciepes for ElectricFields, Simulations, etc.
"""
module Plotting

using ..Common
using RecipesBase
using LaTeXStrings

include("1D/Plotting1D.jl")
# include("2D/Plotting2D.jl")
# include("3D/Plotting3D.jl")

@recipe f(e::ElectricField;by=abs2) = e,by
@recipe f(by::Function,e::ElectricField) = e,by
@recipe f(sim::Simulation,e::ElectricField;by=abs2) = sim,e,by
@recipe f(e::ElectricField,sim::Simulation;by=abs2) = sim,e,by

end
