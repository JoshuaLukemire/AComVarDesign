module AComVarDesign

using Statistics
using LinearAlgebra
using StatsBase
using Test

include("src/Design.jl")
include("src/design_navigation_functions.jl")
include("src/optGA.jl")

export Design, optGA

end
