# MendelKinship

This [Julia](http://julialang.org/) package computes genetic kinship and other identity coefficients. MendelKinship is one component of the umbrella [OpenMendel](https://openmendel.github.io) project.

[![](https://img.shields.io/badge/docs-current-blue.svg)](https://OpenMendel.github.io/MendelKinship.jl)

## Installation

`MendelKinship` currently supports Julia version 1.0 and 1.1, but it is currently an unregistered package. To install, press `]` to invoke the package manager mode and install these packages by typing:
```
add https://github.com/OpenMendel/SnpArrays.jl
add https://github.com/OpenMendel/MendelSearch.jl
add https://github.com/OpenMendel/MendelBase.jl
add https://github.com/biona001/MendelKinship.jl
```

You will also need a few registered packages. Add them by typing:
```
add PlotlyJS Statistics StatsBase CSV ORCA DataFrames
```

## Tutorials
For a tutorial of `MendelKinship`, visit [tutorials](https://github.com/OpenMendel/Tutorials/blob/master/Kinship/KinshipTutorial.ipynb) of the OpenMendel github page. 

## Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ## Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

## Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
