__precompile__()

"""
This module implements computation of kinship and other identity coefficients.
The kinship matrix is computed recursively by the standard formula.
The delta7 matrix is derived from the kinship matrix for non-inbred pedigrees.
The remaining identity coefficients of Jacquard are computed by simulation.
"""
module MendelKinship
#
# Required OpenMendel packages and modules.
#
using MendelBase
using SnpArrays
#
# Required external modules.
#
using CSV
using DataFrames
using ORCA
using PlotlyJS
using StatsBase

export Kinship

"""
This is the wrapper function for the Kinship analysis option.
"""
function Kinship(control_file = ""; args...)

  KINSHIP_VERSION :: VersionNumber = v"0.8.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println("     Kinship analysis option")
  println("        version ", KINSHIP_VERSION)
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  keyword["kinship_output_file"] = ""
  keyword["repetitions"] = 1
  keyword["xlinked_analysis"] = false
  keyword["compare_kinships"] = false
  keyword["maf_threshold"] = 0.01
  keyword["grm_method"] = :MoM # more robust & less sensitive to rare SNPs;
                               # alternatives are :GRM and :Robust
  keyword["deviant_pairs"] = 0
  keyword["kinship_plot"] = ""
  keyword["z_score_plot"] = ""
  keyword["full_pedigree_file"] = ""
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  if keyword["kinship_output_file"] == ""
    keyword["kinship_output_file"] = keyword["output_file"]
  end
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "kinship")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "Kinship"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Check if SNP data were read and names are unique.
  #
  if keyword["compare_kinships"]
    if snpdata.snps == 0
      throw(ArgumentError(
        "\n\nERROR: This analysis requires data from SNP files!\n"))
    elseif snpdata.snps != length(unique(snpdata.snpid))
      throw(ArgumentError(
        " \n\nERROR: There are non-unique names of SNPs.\n"))
    elseif person.people != length(unique(person.name))
      throw(ArgumentError(
        " \n\nERROR: There are non-unique names of individuals.\n"))
    end
  elseif snpdata.snps != 0
      throw(ArgumentError(
        " \n\nERROR: This analysis does not use data from SNP files!\n"))
  end
  #
  # Execute the specified analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = false
  if keyword["compare_kinships"]
    kinship_frame = compare_kinships(pedigree_frame, pedigree, person,
      snpdata, keyword)
  else
    kinship_frame = theoretical_coefficients(pedigree, person, keyword)
  end
  display(kinship_frame)

  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
  end
  #
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing
end # function Kinship

"""
This function orchestrates the comparison of theoretical
and empirical kinships including Fisher z-scores.
The results are placed in a dataframe.
"""
function compare_kinships(pedigree_frame::DataFrame,
  genotyped_pedigree::Pedigree, genotyped_person::Person,
  snpdata::SnpDataStruct, keyword::Dict{AbstractString, Any})
  #
  # Initialize constants
  #
  genotyped_people = genotyped_person.people
  deviant_pairs = keyword["deviant_pairs"]
  if deviant_pairs == 0
    deviant_pairs = genotyped_people^2
  end
  #
  # Record the order of the individuals in the user data upon entry,
  # before MendelBase permutes them.
  # 
  entryorder = convert(Vector{Int64}, pedigree_frame[:EntryOrder])
  #
  # If the full pedigree contains people without SNP information,
  # initialize the full pedigree and person including people without SNPs.
  # The boolean vector `contains_snp` (which is not used)
  # records if there are SNPs: contains_snp[i] means person i has SNP data.
  #
  if keyword["full_pedigree_file"] != ""
    full_pedigree, full_person, contains_snp =
      compute_full_person_and_pedigree(keyword, genotyped_person)
  else
    full_pedigree = genotyped_pedigree
    full_person = genotyped_person
  end
  #
  # Initialize more constants.
  #
  full_pedigrees = full_pedigree.pedigrees
  xlinked = keyword["xlinked_analysis"]
  maf_threshold = keyword["maf_threshold"]
  grm_method = Symbol(keyword["grm_method"])
  #
  # People are reordered, and some may not have SNP info.
  # Compute identification maps to match positions of IDs within two vectors.
  # Recall that snpdata.personid has been permuted; thus, we inverse permute.
  #
  personid_in_snpdata = copy(snpdata.personid)
  invpermute!(personid_in_snpdata, entryorder)
  (name_to_id, id_to_name) = correspond(full_person.name, personid_in_snpdata)
  #
  # Compute the genetic relationship matrix and apply Fisher's transform.
  #
  GRM = grm(snpdata.snpmatrix, method=grm_method, minmaf=maf_threshold)
  clamp!(GRM, -0.99999, 0.99999) 
  GRM .= atanh.(GRM)
  #
  # For each pedigrees, compute empiric kinship coefficients
  # and subtract theoretical kinship coefficients.
  #
  kinship = Vector{Matrix{Float64}}(undef, full_pedigrees)
  for ped = 1:full_pedigrees
    kinship[ped] = kinship_matrix(full_pedigree, full_person, ped, xlinked)
    q = full_pedigree.start[ped] - 1
    for j = 1:full_pedigree.individuals[ped]
      jj = name_to_id[j + q]
      jj == 0 && continue
      for i = j:full_pedigree.individuals[ped]
        ii = name_to_id[i + q]
        ii == 0 && continue
        GRM[ii, jj] = GRM[ii, jj] - atanh(kinship[ped][i, j])
        GRM[jj, ii] = GRM[ii, jj]
      end
    end
  end
  #
  # Find the indices of the most deviant pairs.
  #
  p = partialsortperm(vec(GRM), 1:deviant_pairs, by = abs, rev = true)
  genotyped_indices = CartesianIndices((genotyped_people, genotyped_people))[p]
  #
  # Enter the most deviant pairs in a data frame.
  #
  kinship_frame = DataFrame(Pedigree1 = String[], Pedigree2 = String[],
    Person1 = String[], Person2 = String[], theoretical_kinship= Float64[],
    empiric_kinship=Float64[])
  r = 0.0
  for k = 1:deviant_pairs
    #
    # Find the index of the kth largest deviation in the GRM matrix.
    #
    (ii, jj) = (genotyped_indices[k][1], genotyped_indices[k][2])
    #
    # Map this snpid index back to the individual's name.
    #
    (i, j) = (id_to_name[ii], id_to_name[jj])
    j > i && continue
    (pedi, pedj) = (full_person.pedigree[i], full_person.pedigree[j])
    if pedi == pedj
      q = full_pedigree.start[pedi] - 1
      r = GRM[ii, jj] + atanh(kinship[pedi][i - q, j - q]) 
      push!(kinship_frame,
        [full_pedigree.name[pedi], full_pedigree.name[pedj],
         full_person.name[i], full_person.name[j],
         kinship[pedi][i - q, j - q], tanh(r)])
    else
      push!(kinship_frame,
        [full_pedigree.name[pedi], full_pedigree.name[pedj],
         full_person.name[i], full_person.name[j],
         0.0, tanh(GRM[ii, jj])])
    end
  end
  #
  # Compute the Fisher z-scores.
  #
  kinship_frame = compute_fishers_z(kinship_frame)
  #
  # Create and save plots based on user's requests.
  #
  plot_pair = make_plot_pairs(kinship_frame) 
  if keyword["kinship_plot"] != "" 
    my_compare_plot = make_compare_plot(kinship_frame, plot_pair)
    # PlotlyBase.savefig(my_compare_plot, keyword["kinship_plot"] * ".pdf")
    PlotlyJS.savehtml(my_compare_plot, keyword["kinship_plot"] * ".html")
    println("Kinship plot saved.")
  end
  if keyword["z_score_plot"] != ""
    my_fisher_plot = plot_fisher_z(kinship_frame, plot_pair)
    # PlotlyBase.savefig(my_fisher_plot, keyword["z_score_plot"] * ".pdf")
    PlotlyJS.savehtml(my_fisher_plot, keyword["z_score_plot"] * ".html")
    println("Fisher's plot saved.")
  end
  #
  # Write the combined theoretical and empirical kinship frame to a file.
  #
  kinship_table_file = string(keyword["kinship_output_file"])
  CSV.write(kinship_table_file, kinship_frame;
    writeheader = true, delim = keyword["output_field_separator"],
    missingstring = keyword["output_missing_value"])

  return kinship_frame
end # function compare_kinships

"""
This function orchestrates the computation via gene dropping
of the Kinship and Delta7 coefficients deterministically
and Jacquard's 9 identity coefficients stochastically.
The results are placed in a dataframe.
"""
function theoretical_coefficients(pedigree::Pedigree, person::Person,
  keyword::Dict{AbstractString, Any})

  pedigrees = pedigree.pedigrees
  repetitions = keyword["repetitions"]
  xlinked = keyword["xlinked_analysis"]
  #
  # Define dataframes for theoretical Kinship, Delta7, and Jacquard 
  # coefficients.
  #
  kinship_dataframe = DataFrame(ped1 = Int[], Pedigree = String[],
    Person1 = String[], Person2 = String[], Kinship = Float64[])
  if !xlinked
    delta_dataframe = DataFrame(ped2 = Int[], Pedigree = String[],
      Person1 = String[], Person2 = String[], Delta7 = Float64[])
  end
  coefficient_dataframe = DataFrame(ped3 = Int[], Pedigree = String[],
    Person1 = String[], Person2 = String[],
    delta1 = Float64[], delta2 = Float64[], delta3 = Float64[],
    delta4 = Float64[], delta5 = Float64[], delta6 = Float64[],
    delta7 = Float64[], delta8 = Float64[], delta9 = Float64[])
  #
  # Loop over all pedigrees.
  #
  for ped = 1:pedigrees
    #
    # Retrieve the various coefficients.
    #
    kinship = kinship_matrix(pedigree, person, ped, xlinked)
    if !xlinked
      delta7 = delta7_matrix(pedigree, person, kinship, ped)
    end
    coefficient = jacquard_coefficients(pedigree, person, ped,
      repetitions, xlinked)
    #
    # Insert the computed coefficients into the appropriate dataframes.
    #
    q = pedigree.start[ped] - 1
    for i = 1:pedigree.individuals[ped]
      for j = i:pedigree.individuals[ped]
        push!(kinship_dataframe, [ped, pedigree.name[ped],
          person.name[i + q], person.name[j + q], kinship[i, j]])
        if !xlinked
          push!(delta_dataframe, [ped, pedigree.name[ped],
            person.name[i + q], person.name[j + q], delta7[i, j]])
        end
        push!(coefficient_dataframe, [ped, pedigree.name[ped],
          person.name[i + q], person.name[j + q],
          coefficient[1, i, j], coefficient[2, i, j], coefficient[3, i, j],
          coefficient[4, i, j], coefficient[5, i, j], coefficient[6, i, j],
          coefficient[7, i, j], coefficient[8, i, j], coefficient[9, i, j]])
      end
    end
  end
  #
  # Join the separate dataframes, and sort the combined dataframe so
  # that relative pairs are lumped by pedigree..
  #
  if xlinked
    combined_dataframe = join(kinship_dataframe, coefficient_dataframe,
      on = [:Pedigree, :Person1, :Person2])
    sort!(combined_dataframe, [:ped1, :Person1, :Person2])
    deletecols!(combined_dataframe, [:ped1, :ped3])
  else
    temp_dataframe = join(kinship_dataframe, delta_dataframe,
      on = [:Pedigree, :Person1, :Person2])
    combined_dataframe = join(temp_dataframe, coefficient_dataframe,
      on = [:Pedigree, :Person1, :Person2])
    sort!(combined_dataframe, [:ped1, :Person1, :Person2])
    deletecols!(combined_dataframe, [:ped1, :ped2, :ped3])
  end
  #
  # Write the combined coefficient frame to a file and return.
  #
  kinship_table_file = string(keyword["kinship_output_file"])
  CSV.write(kinship_table_file, combined_dataframe;
    writeheader = true, delim = keyword["output_field_separator"],
    missingstring = keyword["output_missing_value"])
  return combined_dataframe
end # function theoretical_coefficients

"""
This function recursively computes a kinship matrix for each pedigree.
Children must appear after parents. Each pedigree occurs in a block
demarked by its start and finish and is defined implicitly by
the provided fathers and mothers. Person i is a pedigree founder
if and only if mother[i] = 0. X-linked inheritance is allowed.
"""
function kinship_matrix(pedigree::Pedigree, person::Person,
  ped::Int, xlinked::Bool)
  #
  # Allocate a kinship matrix and an offset.
  #
  start = pedigree.start[ped]
  finish = pedigree.finish[ped]
  kinship = zeros(Float64, pedigree.individuals[ped], pedigree.individuals[ped])
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
  cotwin_extension!(kinship, pedigree, person, ped)
  return kinship
end # function kinship_matrix

"""
This function computes the Delta7 matrix for non-inbred pedigrees.
Children must appear after parents. Each pedigree occurs in a block
and is defined implicitly by the provided fathers and mothers.
Person i is a pedigree founder if and only if mother[i] = 0.
"""
function delta7_matrix(pedigree::Pedigree, person::Person,
  kinship::Matrix{Float64}, ped::Int)

  start = pedigree.start[ped]
  finish = pedigree.finish[ped]
  delta7 = zeros(Float64, pedigree.individuals[ped], pedigree.individuals[ped])
  q = start - 1
  for i = start:finish
    for j = start:i
      if j == i
        delta7[i - q, i - q] = 1.0
      elseif person.mother[i] > 0 && person.mother[j] > 0
        k = person.father[i]
        l = person.mother[i]
        m = person.father[j]
        n = person.mother[j]
        delta7[i - q, j - q] = kinship[k - q, m - q] * kinship[l - q, n - q] +
          kinship[k - q, n - q] * kinship[l - q, m - q]
        delta7[j - q, i - q] = delta7[i - q, j - q]
      end
    end
  end
  #
  # Extend the delta7 matrix to co-twins.
  #
  cotwin_extension!(delta7, pedigree, person, ped)
  return delta7
end # function delta7_matrix

"""
This function computes Jacquard's nine condensed identity
coefficients for each pedigree by simulation. Children
must appear after parents. A pedigree is defined implicitly
by giving fathers and mothers. Person i is a pedigree founder
if and only if mother[i] = 0. The more replicates requested,
the more precision attained. X-linked inheritance is possible.
"""
function jacquard_coefficients(pedigree::Pedigree, person::Person,
  ped::Int, replicates::Int, xlinked::Bool)

  start = pedigree.start[ped]
  finish = pedigree.finish[ped]
  ped_size = pedigree.individuals[ped]
  coefficient = zeros(Float64, 9, ped_size, ped_size)
  x = zeros(Float64, ped_size, ped_size)
  source = zeros(Int, 2, ped_size)
  #
  # Loop over all replicates.
  #
  q = start - 1
  for rep = 1:replicates
    f = 0
    for j = start:finish
      if person.mother[j] == 0
        f = f + 1
        source[1, j - q] = f
        if xlinked && person.male[j]
        else
          f = f + 1
        end
        source[2, j - q] = f
      else
        if xlinked && person.male[j]
          u = rand(1:2, 1)
          source[1, j - q] = source[u[1], person.mother[j] - q]
          source[2, j - q] = source[1, j - q]
        else
          u = rand(1:2, 2)
          source[1, j - q] = source[u[1], person.father[j] - q]
          source[2, j - q] = source[u[2], person.mother[j] - q]
        end
      end
    end
    #
    # Loop over all pairs.
    #
    for j = start:finish
      for k = start:finish
        l = identity_state(source[:, j - q], source[:, k - q])
        coefficient[l, j - q, k - q] = coefficient[l, j - q, k - q] + 1.0
      end
    end
  end
  #
  # Extend the coefficients to co-twins.
  #
  for i = 1:9
    x = reshape(coefficient[i, :, :], ped_size, ped_size)
    cotwin_extension!(x, pedigree, person, ped)
    coefficient[i, :, :] = reshape(x, 1, ped_size, ped_size)
  end
  return coefficient / replicates
end # function jacquard_coefficients

"""
This function extends identity coefficients to co-twins.
"""
function cotwin_extension!(x::Matrix{Float64}, pedigree::Pedigree,
  person::Person, ped::Int)

  q = pedigree.start[ped] - 1
  for i = pedigree.finish[ped] + 1:pedigree.twin_finish[ped]
    m = person.primary_twin[i]
    for j = pedigree.start[ped]:pedigree.finish[ped]
      x[i - q, j - q] = x[m - q, j - q]
      x[j - q, i - q] = x[i - q, j - q]
    end
    for j = pedigree.finish[ped] + 1:pedigree.twin_finish[ped]
      n = person.primary_twin[j]
      x[i - q, j - q] = x[m - q, n - q]
    end
  end
end # function cotwin_extension!

"""
This function assigns one of the 9 condensed identity states
to two sourced genotypes. The integer n is the number of
unique genes among the 4 sampled genes.
"""
function identity_state(source1::Vector{Int}, source2::Vector{Int})

  n = length(union(BitSet(source1), BitSet(source2)))
  if n == 4
    return 9
  elseif n == 3
    if source1[1] == source1[2]
      return 4
    elseif source2[1] == source2[2]
      return 6
    else
      return 8
    end
  elseif n == 2
    if source1[1] == source1[2]
      if source2[1] == source2[2]
        return 2
      else
        return 3
      end
    else
      if source2[1] == source2[2]
        return 5
      else
        return 7
      end
    end
  else
    return 1
  end
end # function identity_state


"""
This function assembles pairs of names and values for printing.
"""
function make_plot_pairs(x::DataFrame)

  plot_pair = Array{String}(undef, size(x, 1))
  for i in 1:length(plot_pair)
    plot_pair[i] = "Person1=" * x[i, 3] * ", " * "Person2=" * x[i, 4]
  end
  return plot_pair
end # function make_plot_pairs


"""
This function composes the comparison plot using PlotlyJS.
"""
function make_compare_plot(x::DataFrame, plot_pair::Vector{String})

  trace1 = scatter(;x=x[:theoretical_kinship], 
    y=x[:empiric_kinship], mode="markers", 
    name="empiric kinship", text=plot_pair)
  trace2 = scatter(;x=[1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 0.0],
    y=[1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 0.0], 
    mode="markers", name="marker for midpoint")
  layout = Layout(;title="Compare Empiric and Theoretical Kinships",
    hovermode="closest", 
    xaxis=attr(title="Theoretical Kinship", 
      showgrid=false, zeroline=false),
    yaxis=attr(title="Empiric Kinship", zeroline=false))
  data = [trace1, trace2]
  plot(data, layout)
end # function make_compare_plot


"""
This function computes Fisher's atanh transform and standardizes.
"""
function compute_fishers_z(x::DataFrame)

  theoretical_transformed = map(atanh, x[:theoretical_kinship])
  empiric_transformed = map(atanh, x[:empiric_kinship])
  difference = empiric_transformed - theoretical_transformed
  new_zscore = zscore(difference) #in StatsBase package
  return [x DataFrame(fishers_zscore = new_zscore)]
end # function compute_fishers_z


"""
This function plots the transformed data using PlotlyJS.
"""
function plot_fisher_z(x::DataFrame, plot_pair::Vector{String})

    trace1 = histogram(x=x[:fishers_zscore], text=plot_pair)
    data = [trace1]
    layout = Layout(barmode="overlay", 
        title="Z-score plot for Fisher's statistic",
        xaxis=attr(title="Standard Deviations"),
        yaxis=attr(title="Count"))
    plot(data, layout)
end # function plot_fisher_z


"""
 
"""
function correspond(x::AbstractVector{AbstractString},
  y::AbstractVector{AbstractString})

  (m, n) = (length(x), length(y))
  xperm = sortperm(x)
  yperm = sortperm(y)
  x_to_y = zeros(Int, m)
  y_to_x = zeros(Int, n)
  (i, j) = (1, 1)
  done = false
  while !done
    ii = xperm[i]
    jj = yperm[j]
    if x[ii] == y[jj]
      x_to_y[ii] = jj
      y_to_x[jj] = ii
      (i, j) = (i + 1, j + 1)
    elseif x[ii] < y[jj]
      i = i +1
    elseif y[jj] < x[ii]
      j = j + 1
    end
    done = i > m || j > n
  end
  return (x_to_y, y_to_x)
end # function correspond

# x = ["the", "and", "a", "be", "an"]
# y = ["but", "be", "a", "an", "so"]
# (x_to_y, y_to_x) = correspond(x, y)


"""

"""
function compute_full_person_and_pedigree(keyword :: Dict{AbstractString, Any},
  subset_person :: Person)
#
# Use MendelBase's function to construct the "full" pedigree and person data.
# First, store original key values.
#
  original_pedigree_file = keyword["pedigree_file"]
  snpdefinition_file = keyword["snpdefinition_file"]
  snpdata_file = keyword["snpdata_file"]
  
  keyword["pedigree_file"] = keyword["full_pedigree_file"]
  keyword["snpdefinition_file"] = "" # Avoid reading in snpdata again.
  keyword["snpdata_file"] = ""
  full_pedigree, full_person, = read_external_data_files(keyword)
  keyword["pedigree_file"] = original_pedigree_file
  keyword["snpdefinition_file"] = snpdefinition_file
  keyword["snpdata_file"] = snpdata_file
#
# Now compute indicator of person i having SNP information.
#
  total_people = full_person.people
  contains_snp = falses(total_people)
  for i in 1:total_people
    person_i = full_person.name[i]
    (person_i in subset_person.name) && (contains_snp[i] = true)
  end

  return full_pedigree, full_person, contains_snp
end # function compute_full_person_and_pedigree

#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module MendelKinship
