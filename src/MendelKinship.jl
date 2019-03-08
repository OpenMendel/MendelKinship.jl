module MendelKinship

using SnpArrays
using MendelBase
using PlotlyJS
using StatsBase
using DataFrames
using CSV
using ORCA

export Kinship
export compare_kinships, kinship_option, kinship_matrix, delta7_matrix, jacquard_coefficients
export cotwin_extension!, identity_state, correspond, compute_full_pedigree

include("kinship_utilities.jl")

end # module MendelKinship
