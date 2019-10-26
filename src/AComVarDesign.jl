module AComVarDesign

using Statistics
using LinearAlgebra
using StatsBase
using Test

include("Design.jl")
include("design_navigation_functions.jl")
include("optGA.jl")

export Design, optGA

end
