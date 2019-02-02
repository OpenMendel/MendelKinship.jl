module MendelKinship

# using MendelBase
using SnpArrays
using PlotlyJS
using StatsBase
using DataFrames
using CSV

# MendelBase dependencies which would be removed after MendelBase is at 1.0 on github
using Distributions
using GLM
using LinearAlgebra
using Missings
using Printf
using Random
using Statistics

export Kinship

include("MendelBase.jl")
include("kinship_utilities.jl")

end # module MendelKinship

