*Julia utilities for computing theoretical/empirical kinship and other identity coefficients*

Mendel Kinship is a component of the umbrella [OpenMendel project](https://github.com/OpenMendel). 

## Package features
+ Calculation of empirical and theoretical kinship coefficients, and their comparisons
+ Calculation of Jacquard's 9 condensed identity coefficients
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
    "man/api.md",
]
Depth = 2
```
