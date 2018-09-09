using MendelKinship
using MendelBase
using DataFrames 

srand(123)

@testset "kinship_matrix" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["kinship_output_file"] = "Kinship_Output_File.txt"
    keyword["compare_kinships"] = false
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
        b = mod(b - 1, 4) + 1 # subtract 1, mod, then add one. Doing so skips 0.
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
    keyword["compare_kinships"] = false
    process_keywords!(keyword, "kinship_Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)
    #
    #   This is the probability that two individuals share both of their alleles 
    #   identity by descent. It's case 7 on page 10 of 
    #   http://faculty.washington.edu/tathornt/BIOST551/lectures_2012/Lecture7_Coefficients_of_Identity_and_Coancestry.pdf
    #
    kinship_matrix = MendelKinship.kinship_matrix(pedigree, person, 1, false)
    delta7_matrix = MendelKinship.delta7_matrix(pedigree, person, kinship_matrix, 1)
    @test issymmetric(delta7_matrix)
    @test eltype(delta7_matrix) == Float64
    @test size(delta7_matrix) == (6, 6)

    for i = 1:size(delta7_matrix, 1)
        @test delta7_matrix[i, i] == 1.0 #one is always IBD for both alleles to oneself
    end

    # founders are not related, and founders and their kid cannot be fraternity brothers
    # since 1 parent can only pass on 1 of his alleles. 
    @test delta7_matrix[1, 2] == 0.0 
    @test delta7_matrix[1, 3] == 0.0
    @test delta7_matrix[1, 4] == 0.0
    @test delta7_matrix[1, 5] == 0.0
    @test delta7_matrix[1, 6] == 0.0

    # Assume father = (A1, A2) and mother = (A3, A4), with 2 kids = person 3 and 4.
    # Say father passes down A1 to person 3, then person 4 also receiving this has 
    # probability 0.5. Then the mother passes down A3 or A4, so the other kid also receiving 
    # it has probability 0.5. Thus person 3 & 4 shares alleles with probability 0.5 * 0.5. 
    @test delta7_matrix[3, 4] == 0.5 * 0.5
    # Then person 3 and person 4 mated, giving birth to person 5 & 6. Person 3 donates one of 
    # its allele with probability 1, and for person 5 to be share both alleles with its parents,
    # person 4 must donate the correct half with probability 0.5. Thus
    @test delta7_matrix[3, 5] == 0.5 * (0.5 * 0.5)
    @test delta7_matrix[3, 6] == 0.5 * (0.5 * 0.5)
    @test delta7_matrix[4, 5] == 0.5 * (0.5 * 0.5)
    @test delta7_matrix[4, 6] == 0.5 * (0.5 * 0.5)
    # Finally,
    # P(5 & 6 share 2 alleles) = 
    #     P(5 and 6 share 2 | 3 & 4 share none) 
    #   + P(5 and 6 share 2 | 3 & 4 share 1) 
    #   + P(5 and 6 share 2 | 3 & 4 share 2) ... which we can calculate by hand:
    @test delta7_matrix[5, 6] == 0.25 * 0.25 + 0.5 * 0.25 + 0.25 * 0.5
end

#
# Need a function that spits out a kinship matrix without extening to co-twins
# to fully test whether cotwin_extension is working properly.
#
function kinship_matrix_notwin(pedigree::Pedigree, person::Person,
  ped::Int, xlinked::Bool)
  #
  # Allocate a kinship matrix and an offset.
  #
  start = pedigree.start[ped]
  finish = pedigree.finish[ped]
  kinship = zeros(pedigree.individuals[ped], pedigree.individuals[ped])
  q = start - 1
  for i = start:finish
    #
    # Initialize the kinship coefficients for founders.
    #
    j = person.mother[i]
    if j == 0
      if xlinked && person.male[i]
        kinship[i - q, i - q] = 1.0
      else
        kinship[i - q, i - q] = 0.5
      end
    else
      #
      # Compute the kinship coefficients of descendants by averaging.
      #
      if xlinked && person.male[i]
        kinship[i - q, i - q] = 1.0
        for m = start:i - 1
          kinship[i - q, m - q] = kinship[j - q, m - q]
          kinship[m - q, i - q] = kinship[i - q, m - q]
        end
      else
        k = person.father[i]
        kinship[i - q, i - q] = 0.5 + 0.5 * kinship[j - q, k - q]
        for m = start:i - 1
          kinship[i - q, m - q] = 0.5 * (kinship[j - q, m - q]
             + kinship[k - q, m - q])
          kinship[m - q, i - q] = kinship[i - q, m - q]
        end
      end
    end
  end
  #
  # Extend the kinship matrix to co-twins.
  #
  # cotwin_extension!(kinship, pedigree, person, ped)
  return kinship
end # function kinship_matrix_notwin

@testset "cotwin_extension!" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["kinship_output_file"] = "Kinship_Output_File.txt"
    keyword["compare_kinships"] = false
    process_keywords!(keyword, "kinship_Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    # compute kinship matrix without extending to cotwins
    matrix1 = kinship_matrix_notwin(pedigree, person, 1, false) 
    matrix2 = kinship_matrix_notwin(pedigree, person, 2, false)

    @test issymmetric(matrix1)
    @test issymmetric(matrix2)
    @test eltype(matrix2) == Float64
    @test eltype(matrix1) == Float64
    @test size(matrix1) == (6, 6) 
    @test size(matrix2) == (7, 7)

    @test length(find(matrix1)) == 34 # only 2 element of 6x6 matrix is 0
    @test all(matrix2[6, :] .== 0.0) == true #twins have 0 values before calling cotwin_extension
    @test all(matrix2[7, :] .== 0.0) == true

    MendelKinship.cotwin_extension!(matrix1, pedigree, person, 1)
    MendelKinship.cotwin_extension!(matrix2, pedigree, person, 2)

    @test size(matrix1) == (6, 6) # size didn't change
    @test size(matrix2) == (7, 7)
    @test length(find(matrix1)) == 34 # non-zero values should not change since no twins
    @test matrix2[6, :] == matrix2[5, :] #twins got reassigned
    @test matrix2[7, :] == matrix2[5, :] #twins got reassigned
end

@testset "jacquard_coefficients" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["kinship_output_file"] = "Kinship_Output_File.txt"
    keyword["compare_kinships"] = false
    process_keywords!(keyword, "kinship_Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    n = 1000 #default number of repetition
    matrix = MendelKinship.jacquard_coefficients(pedigree, person, 1, n, false) 
    @test size(matrix) == (9, 6, 6)
    @test eltype(matrix) == Float64

    # testing some values, but they weren't checked by hand.
    @test all(matrix[1, :, 1] .== 0.0)
    @test matrix[3, 5, 1] == 0.119
    @test matrix[8, 6, 1] == 0.491
    @test all(matrix[2, :, 2] .== 0.0)
    @test matrix[4, 5, 2] == 0.119
    @test matrix[9, 1, 2] == 1.0
    @test matrix[3, 3, 3] == 0.0
    @test matrix[7, 4, 3] == 0.23
    @test matrix[8, 5, 3] == 0.509
    @test matrix[3, 6, 4] == 0.229
    @test all(matrix[5, :, 4] .== 0.0)
    @test matrix[1, 5, 5] == 0.237
    @test matrix[1, 4, 5] == 0.0
    @test matrix[5, 1, 5] == 0.119
    @test matrix[6, 6, 6] == 0.0
    @test matrix[2, 5, 6] == 0.031
    @test matrix[7, 2, 6] == 0.128
end

@testset "basics & wrapper functions" begin
    keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
    keyword["kinship_output_file"] = "Kinship_Output_File.txt"
    keyword["xlinked_analysis"] = false
    keyword["compare_kinships"] = false
    process_keywords!(keyword, "kinship_Control.txt", "")
    (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
        read_external_data_files(keyword)

    result = MendelKinship.kinship_option(pedigree, person, keyword)
    kinship_matrix = MendelKinship.kinship_matrix(pedigree, person, 1, false)
    delta7_matrix = MendelKinship.delta7_matrix(pedigree, person, kinship_matrix, 1)

    # test values to see if numbers are matching. 
    # Note that result matrix removed all redundant rows,
    # (e.g. X[3, 4] == X[4, 3] so one X[3, 4] will be listed)
    @test size(result) == (49, 14)
    @test all(result[1:6, 5] .== delta7_matrix[1, 1:6])
    @test all(result[7:11, 5] .== delta7_matrix[2, 2:6]) 
    # @test all(result[12:15, 5] .== delta7_matrix[3, 3:6]) #turns out the structure of output is more compact than this.
    # @test all(result[16:18, 5] .== delta7_matrix[4, 4:6]) 
    # @test all(result[19:20, 5] .== delta7_matrix[5, 5:6]) 
    #@test result[21, 1] == delta7_matrix[6, 6]
    @test all(result[1:6, 4] .== kinship_matrix[1, 1:6])
    @test all(result[7:11, 4] .== kinship_matrix[2, 2:6]) 
    #@test all(result[12:15, 4] .== kinship_matrix[3, 3:6]) #turns out the structure of output is more compact than this.
    #@test all(result[16:18, 4] .== kinship_matrix[4, 4:6]) 
    #@test all(result[19:20, 4] .== kinship_matrix[5, 5:6]) 
    #@test result[21, 1] == kinship_matrix[6, 6]

    #final test, if no error then should return nothing.
    @test Kinship("kinship_Control.txt") == nothing

    #current coverage = 145/158
end