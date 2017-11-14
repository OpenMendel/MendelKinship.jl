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
    # 2. Φii = 0.5 + 0.5Φkl when i not founder, and k, l are parents of i
    # 3. Φij = 0 if i, j both founders (i.e. founders are not related)
    # 4. Φij = 0.5 Φjk + 0.5Φjl when i not founder (i>j), k and l are parents of i
    # 5. Parents precedes children, and either both parent are present or both aren't, but 
    #    this was tested in MnedelBase. 
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
end

@testset "identity_state" begin
    #
    # Testing all possible combinations (hence the weird modular arithmetics). 
    # Basically match up geometrically on page 10:
    # http://faculty.washington.edu/tathornt/BIOST551/lectures_2012/Lecture7_Coefficients_of_Identity_and_Coancestry.pdf
    #
    for i = 1:4
        a, b, c, d = i, i + 1, i + 2, i + 3
        b = mod(b - 1, 4) + 1 # subtract 1, mod, then add one at the end. This skips 0.
        c = mod(c - 1, 4) + 1
        d = mod(d - 1, 4) + 1

        @test MendelKinship.identity_state([a, a], [a, a]) == 1
        @test MendelKinship.identity_state([a, a], [b, b]) == 2
        @test MendelKinship.identity_state([a, a], [a, b]) == 3
        @test MendelKinship.identity_state([a, a], [b, a]) == 3
        @test MendelKinship.identity_state([a, a], [b, c]) == 4
        @test MendelKinship.identity_state([a, b], [a, a]) == 5
        @test MendelKinship.identity_state([a, b], [b, b]) == 5
        @test MendelKinship.identity_state([a, b], [c, c]) == 6
        @test MendelKinship.identity_state([a, b], [a, b]) == 7
        @test MendelKinship.identity_state([a, b], [b, a]) == 7
        @test MendelKinship.identity_state([a, b], [a, c]) == 8
        @test MendelKinship.identity_state([a, b], [c, a]) == 8
        @test MendelKinship.identity_state([a, b], [b, c]) == 8
        @test MendelKinship.identity_state([a, b], [c, b]) == 8
        @test MendelKinship.identity_state([a, b], [c, d]) == 9
    end
end

@testset "delta7_matrix" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["kinship_output_file"] = "Kinship_Output_File.txt"
    process_keywords!(keyword, "kinship_Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    kinship_matrix = MendelKinship.kinship_matrix(pedigree, person, 1, false)
    delta7_matrix = MendelKinship.delta7_matrix(pedigree, person, kinship_matrix, 1)
    #
    #   This seems to be the fraternity coefficient. It is the probability that 
    #   two individuals share both of their alleles identity by descent.
    #


end

@testset "jacquard_coefficients" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["kinship_output_file"] = "Kinship_Output_File.txt"
    process_keywords!(keyword, "kinship_Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    MendelKinship.jacquard_coefficients(pedigree, person, 1, 1000, false) #1000 repetition = default


end

@testset "cotwin_extension!" begin
    
end

@testset "basics & wrapper functions" begin
    
end