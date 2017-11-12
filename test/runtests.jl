module PkgTest
using MendelKinship
using Base.Test

# write your own tests here
include("MendelKinship_test.jl")


# julia -e 'Pkg.test("MendelKinship",coverage=true)'
# @show get_summary(process_file("src/MendelKinship.jl"))

