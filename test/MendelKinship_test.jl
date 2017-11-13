using MendelKinship
using MendelBase
using DataFrames 

@testset "kinship_matrix" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["kinship_output_file"] = "Kinship_Output_File.txt"
    process_keywords!(keyword, "kinship_Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    @test_throws(MethodError, MendelKinship.kinship_matrix("random", person, 1, false))
    @test_throws(MethodError, MendelKinship.kinship_matrix(pedigree, "person", 1, false))
    @test_throws(BoundsError, MendelKinship.kinship_matrix(pedigree, person, 3, false)) #only 2 pedigree
    @test_throws(MethodError, MendelKinship.kinship_matrix(pedigree, person, Inf, false))
    @test_throws(MethodError, MendelKinship.kinship_matrix(pedigree, person, NaN, false))
    @test_throws(BoundsError, MendelKinship.kinship_matrix(pedigree, person, 0, false))    

    #
    # The following rules define the kinship coefficient matrix:
    # 1. Φii = 0.5 if i is a founder
    # 2. Φij = 0 if i, j both founders (i.e. founders are not related)
    # 3. Φij = 0.5 Φjk + 0.5Φjl when i not founder (i>j), k and l are parents of i
    # 4. Φii = 0.5 + 0.5Φkl when i not founder, and k, l are parents of i
    #
    matrix = MendelKinship.kinship_matrix(pedigree, person, 1, false)

    @test size(matrix) == (6, 6)
    @test eltype(matrix) == Float64
    @test issymmetric(matrix) #kinship matrices are always symmetric

    #founders
    @test matrix[1, 1] == 0.5 
    @test matrix[2, 2] == 0.5 
    @test matrix[1, 2] == 0.0 

    #first row & first column
    @test matrix[3, 1] == 0.5 * matrix[1, 2] + 0.5 * matrix[1, 1] 
    @test matrix[4, 1] == 0.5 * matrix[1, 2] + 0.5 * matrix[1, 1] 
    @test matrix[5, 1] == 0.5 * matrix[1, 3] + 0.5 * matrix[1, 4]
    @test matrix[6, 1] == 0.5 * matrix[1, 3] + 0.5 * matrix[1, 4]

    #second row & second column
    @test matrix[3, 2] == 0.5 * matrix[2, 1] + 0.5 * matrix[2, 2]
    @test matrix[4, 2] == 0.5 * matrix[2, 1] + 0.5 * matrix[2, 2]
    @test matrix[5, 2] == 0.5 * matrix[2, 3] + 0.5 * matrix[2, 4]
    @test matrix[6, 2] == 0.5 * matrix[2, 3] + 0.5 * matrix[2, 4]

    #thrid row & third column
    @test matrix[3, 3] == 0.5 + 0.5 * matrix[1, 2] 
    @test matrix[4, 3] == 0.5 * matrix[3, 1] + 0.5 * matrix[3, 2]
    @test matrix[5, 3] == 0.5 * matrix[3, 3] + 0.5 * matrix[3, 4]
    @test matrix[6, 3] == 0.5 * matrix[3, 3] + 0.5 * matrix[3, 4]

    #4th row & 4th column
    @test matrix[4, 4] == 0.5 + 0.5 * matrix[1, 2] 
    @test matrix[5, 4] == 0.5 * matrix[4, 3] + 0.5 * matrix[4, 4]
    @test matrix[6, 4] == 0.5 * matrix[4, 3] + 0.5 * matrix[4, 4]

    #5th row & 5th column
    @test matrix[5, 5] == 0.5 + 0.5 * matrix[3, 4] 
    @test matrix[6, 5] == 0.5 * matrix[5, 3] + 0.5 * matrix[5, 4]

    #6th row & 6th column
    @test matrix[6, 6] == 0.5 + 0.5 * matrix[3, 4] 

    #Too lazy to write tests for another pedigree
    #matrix2 = MendelKinship.kinship_matrix(pedigree, person, 2, false)
end

@testset "delta7_matrix" begin
    
end

@testset "jacquard_coefficients" begin
    
end

@testset "cotwin_extension!" begin
    
end

@testset "identity_state" begin
    
end

@testset "basics & wrapper functions" begin
    
end