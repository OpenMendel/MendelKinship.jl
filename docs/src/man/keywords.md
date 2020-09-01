
# Analysis keywords

| Keyword | Default Value | Allowed value | Description |
| --- | --- | --- | --- |
|`kinship_output_file` | Kinship_Output_File.txt | true/false | OpenMendel generated output file with table of kinship coefficients |
|`repetitions` | 1 | Integer | Repetitions for sharing statistics |
|`xlinked_analysis` | false | beelean| Whether markers are on the X chromosome |
|`compare_kinships` | false | boolean | Whether we want to compare theoretical vs empiric kinship |
|`kinship_plot` | "" | User defined file name | A user specified name for a plot comparing theoretical and empiric kinship value |
|`z_score_plot` | "" | User defined file name |  A user specified name for a plot of fisher's z statistic.  |
|`grm_method` | MoM | GRM, MoM, Robust | Method used for empiric kinship calculation. Defaults to `MoM`, but user could choose the more common `GRM` or Robust GRM methods instead. (**Warning:** Based on our experience, Fisher's z score is very unreliable if the GRM method is used for rare (maf < 0.2) snps) |
|`maf_threshold` | 0.01 | Real number between 0 and 1 | The minor allele frequency threshold for the GRM computation |
|`deviant_pairs` | false | Integer less than $n(n+1)/2$ | Number of top deviant pairs (theoretical vs empiric kinship) the user wants to keep |

## See also

+ [MendelBase documentation](https://openmendel.github.io/MendelBase.jl/#keywords-table) for general keywords common to most analysis package
+ [Old MendelKinship documentation](https://openmendel.github.io/MendelKinship.jl/)
