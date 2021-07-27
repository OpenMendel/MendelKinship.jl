using MendelKinship
using Test
using DataFrames
using CSV

datadir = joinpath(@__DIR__, "../data/documentation_data/")
cd(datadir)

@testset "Example 2 in Tutorial" begin
    Kinship("control_compare_29a.txt")
    @test isfile("kinship_file_output.txt")
 
    df = CSV.read("kinship_file_output.txt", DataFrame) # read result
    @test size(df) == (22578, 7)
    @test df[1, 1] == 14
    @test df[1, 2] == 14
    @test df[1, 3] == 26732
    @test df[1, 4] == 264
    @test df[1, 5] == 0.0
    @test isapprox(df[1, 6], 0.10955177290095937, atol=1e-7)
    @test isapprox(df[1, 7], 5.319415124407314, atol=1e-7)
end
