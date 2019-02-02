################################################################################
# This set of functions provides various genetic utilities.
################################################################################
#
# Required external packages.
#
# using Statistics, Distributions

export map_function, inverse_map_function
export hardy_weinberg_test, xlinked_hardy_weinberg_test

"""
Calculate recombination fractions based on Haldane's or Kosambi's formula.
"""
function map_function(d::T where T <: Real, choice::String)
#
  one_d = one(typeof(d))
  if choice == "Haldane"
    theta = (one_d - exp(-2d)) / 2
  elseif choice == "Kosambi"
    theta = tanh(2d) / 2
  else
  throw(ArgumentError(
    "Only Haldane's or Kosambi's map function is allowed.\n \n"))
  end
  return theta
end # function map_function

"""
Calculates genetic map distances based on Haldane's or Kosambi's formula.
"""
function inverse_map_function(theta::T where T <: Real, choice::String)
#
  one_th = one(typeof(theta))
  if choice == "Haldane"
    if theta >= one_th / 2
      d = Inf
    else
      d = -(log(one_th - 2 * theta)) / 2
    end
  elseif choice == "Kosambi"
    if theta >= one_th / 2
      d = Inf
    else
      d = log((one_th + 2 * theta) / (one_th - 2 * theta)) / 4
    end
  else
    throw(ArgumentError(
      "Only Haldane's or Kosambi's map function is allowed.\n \n"))
  end
  return d
end # function inverse_map_function

"""
Test for Hardy-Weinberg equilibrium at a SNP.
The genotype vector conveys for each person his or her number of
reference alleles. Thus, all genotypes belong to {0, 1, 2, Missing}.
The pvalue of the chi-square statistic is returned.
"""
function hardy_weinberg_test(genotype::Vector)
#
  (expected, observed) = (zeros(3), zeros(3))
  #
  # Tally observations per cell.
  #
  n = 0
  for i = 1:length(genotype)
    if !ismissing(genotype[i])
      n = n + 1
      g = round(Int, genotype[i])
      j = g + 1
      observed[j] = observed[j] + 1.0
    end
  end
  if n == 0
    return 1.0
  end
  #
  # Compute allele frequencies.
  #
  p = (2.0 * observed[1] + observed[2]) / (2 * n)
  if p <= 0.0 || p >= 1.0
    return 1.0
  else
    q = 1.0 - p
    #
    # Sum over cells of the contingency table.
    #
    expected[1] = n * p^2
    expected[2] = n * 2.0 * p * q
    expected[3] = n * q^2
    test = 0.0
    for i = 1:3
    test = test + (observed[i] - expected[i])^2 / expected[i]
  end
  return ccdf(Chisq(1), test)
  end
end # function hardy_weinberg_test

"""
Test for Hardy-Weinberg equilibrium at an xlinked SNP.
The genotype vector conveys for each person his or her number of
reference alleles. Thus, all genotypes belong to {0, 1, 2, Missing}.
The p-value of the chi-square statistic is returned.
"""
function xlinked_hardy_weinberg_test(genotype::Vector, 
  male::BitArray{1})

  (expected_female, observed_female) = (zeros(3), zeros(3))
  (expected_male, observed_male) = (zeros(2), zeros(2))
  #
  # Tally observations per cell and count members of each sex.
  #
  (females, males) = (0, 0)
  for i = 1:length(genotype)
    if !ismissing(genotype[i])
      g = round(Int, genotype[i])
      if male[i]
        males = males + 1
        j = min(g + 1, 2)
        observed_male[j] = observed_male[j] + 1.0
      else
        females = females + 1
        j = g + 1
        observed_female[j] = observed_female[j] + 1.0
      end
    end
  end
  if females + males == 0
    return 1.0
  end
  #
  # Compute allele frequencies.
  #
  p = 2.0 * observed_female[1] + observed_female[2] + observed_male[1]
  p = p / (2 * females + males)
  if p <= 0.0 || p >= 1.0
    return 1.0
  end
  q = 1.0 - p
  #
  # Sum over female cells of the contingency table.
  #
  expected_female[1] = females * p^2
  expected_female[2] = females * 2.0 * p * q
  expected_female[3] = females * q^2
  test = 0.0
  for i = 1:3
    term = (observed_female[i] - expected_female[i])^2 / expected_female[i]
    test = test + term
  end
  #
  # Sum over male cells of the contingency table.
  #
  expected_male[1] = males * p
  expected_male[2] = males * q
  for i = 1:2
    term = (observed_male[i] - expected_male[i])^2 / expected_male[i]
    test = test + term
  end
  return ccdf(Chisq(2), test)
end # function xlinked_hardy_weinberg_test

# Testing Routines
#
# d = rand()
# choice = "Haldane";
# theta = map_function(d, choice)
# di = inverse_map_function(theta, choice);
# println(di);
# d = rand()
# choice = "Kosambi";
# theta = map_function(d, choice)
# di= inverse_map_function(theta, choice);
# println(di);
# p = 1 / 3;
# n = 100;
# x = zeros(n);
# for i = 1:n
#   u = rand(2)
#   if u[1] <= p; x[i] = x[i] + 1.0; end
#   if u[2] <= p; x[i] = x[i] + 1.0; end 
# end
# pvalue = hardy_weinberg_test(x);
# println("pvalue = ",pvalue);
# male = falses(n);
# male[1:div(n, 2)] .= true;
# for i = 1:div(n, 2)
#   x[i] = 0.0
#   u = rand(1)
#   if u[1] <= p; x[i] = x[i] + 1.0; end
# end
# pvalue = xlinked_hardy_weinberg_test(x, male);
# println("pvalue = ",pvalue);
