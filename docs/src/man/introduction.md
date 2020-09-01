
# About

This is a Julia implementation of option 29 of the [Mendel software](https://academic.oup.com/bioinformatics/article/29/12/1568/292810), originally programmed in Fortran. `MendelKinship.jl` contains new features, but its usage largely resembles Mendel by employing *control files* - requiring users to make a file to specify options instead of interacting with Julia directly. This is not the most Julian way of programming or usage experience... but unfortuantely we lack manpower to correct it. 

## When to use MendelKinship

`MendelKinship.jl` is capable of calculating the theoretical kinship coefficient $\Phi_{ij}$ as long as a [valid pedigree structure](https://openmendel.github.io/MendelBase.jl/#pedigree-file) is provided. When SNP markers are available, `MendelKinship.jl` can also calculate empirical kinship coefficients using GRM, robust GRM, or Method of Moments methods (see [this paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20584) and [this paper](https://academic.oup.com/bioinformatics/article/26/22/2867/228512) for details). Here we recommend the Robust GRM or MoM (default) method because their estimates are more robust in the presence of rare alleles. 

`MendelKinship` can optionally compare the empirical kinship and theoretical kinships to check suspect pedigree structures and reveal hidden relatedness. It can also reveal sample mixed ups or other laboratory errors that can lead to inaccurate empirical kinships.  The result is saved in a table sorted in descending order. We optioanlly output 2 interactive plots that allow users to quickly pinpoint pairs with the greatest theoretical vs empirical deviance. 

## Using PLINK files as input

MendelKinship accepts [PLINK binary format](https://www.cog-genomics.org/plink2/formats#bed) as input, in which case the triplets (`data.bim`, `data.bed`, `data.fam`) must all be present. One can import the data by specifying the following in the control file:

`plink_input_basename = data` 

However, sometimes the `.fam` file contains non-unique person id (2nd column of `.fam` file) across different pedigrees, which is currently **not** permitted in MendelKinship. A person's id cannot be repeated in other pedigrees, even if it is contextually clear that they are different persons.
