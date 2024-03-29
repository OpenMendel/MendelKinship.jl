# MendelKinship.jl 

*Julia utilities for computing theoretical/empirical kinship and other identity coefficients*

`MendelKinship.jl` is a component of the umbrella [OpenMendel project](https://link.springer.com/article/10.1007/s00439-019-02001-z). Please read section on **Kinship Comparison**. 

## Package features

+ Calculation of empirical and theoretical kinship coefficients, and their comparisons
+ Approximation of Jacquard's 9 condensed identity coefficients via Monte Carlo gene dropping
+ Calculation of Delta7 matrix for non-inbred pedigrees

For details, please read Chapter 5 of **Mathematical and Statistical Methods for Genetic Analysis** by Kenneth Lange, who wrote most of the code for this package. 

## Installation

Press `]` to invoke the package manager mode and install the following:
```
add https://github.com/OpenMendel/SnpArrays.jl
add https://github.com/OpenMendel/MendelSearch.jl
add https://github.com/OpenMendel/MendelBase.jl
add https://github.com/OpenMendel/MendelKinship.jl
```

This package supports Julia `v1.0`+.

## Manual Outline

```@contents
Pages = [
    "man/introduction.md",
    "man/KinshipTutorial.md"
    "man/keywords.md",
    "man/api.md"
]
Depth = 2
```
