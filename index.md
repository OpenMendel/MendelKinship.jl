### Overview
Mendel Kinship is a component of the umbrella [OpenMendel](https://openmendel.github.io) project. Kinship coefficients quantify the degree of relationship between two relatives. Mendel Kinship computes global kinship coefficients and local kinship coefficients. (Global kinship coefficients are also known as theoretical; local are also known as empirical or conditional.) A global kinship coefficient between two individuals *i* and *j* is the probability that a randomly sampled gene from *i* is identical by descent to a randomly sampled gene from the same arbitrary locus of *j*.

### Appropriate Problems and Data Sets
Global kinship coefficients based solely on pedigree structure can be quickly computed on pedigrees of all sizes, since marker data are ignored. For any local or SNP-based kinship analysis, avoid mixing autosomal and X-linked loci. Both sets of loci can be analyzed, but they must be run separately. For standard, text-based input files, local kinship analyses will be likelihood-based. For such analyses, including comparing the local to global kinship coefficients, large pedigrees will be bypassed. For binary data files, the SNP-based global and local methods will be used.

### Installation
*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [Search](https://openmendel.github.io/Search.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install MendelKinship:

    Pkg.clone("https://github.com/OpenMendel/MendelKinship.jl.git")

This package supports Julia v0.4 and v0.5.

### Input Files
The Mendel Kinship analysis package uses the following input files. Example input files can be found in the [docs]( https://github.com/OpenMendel/MendelKinship.jl/tree/master/docs) subfolder of the Mendel Kinship project. (An analysis won't always need every file type below.)

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](https://openmendel.github.io/MendelBase.jl/#keywords-table)).
* [Locus File](https://openmendel.github.io/MendelBase.jl/#locus-file): Names and describes the genetic loci in your data.
* [Pedigree File](https://openmendel.github.io/MendelBase.jl/#pedigree-file): Gives information about your individuals, such as name, sex, family structure, and ancestry.
* [Phenotype File](https://openmendel.github.io/MendelBase.jl/#phenotype-file): Lists the available phenotypes.
* [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, allele names, allele frequencies.
* [SNP Data File](https://openmendel.github.io/MendelBase.jl/#snp-data-file): Holds the genotypes for your data set. Must be a standard binary PLINK BED file in SNP major format. If you have a SNP data file you must have a SNP definition file.

### Control file<a id="control-file"></a>
The Control file is a text file consisting of keywords and their assigned values. The format of the Control file is:

	Keyword = Keyword_Value(s)

Below is an example of a simple Control file to run Kinship:

	#
	# Input and Output files.
	#
	pedigree_file = kinship_PedigreeFrame.txt
	kinship_output_file = kinship_file_Output.txt
	output_file = kinship_Output.txt
	#
	# Analysis parameters for Kinship option.
	#
	repetitions = 1000

In the example above, there are four keywords. The first keyword specifies the input Pedigree file: *kinship_PedigreeFrame.txt*. The next two keywords specify output files with results of the analysis: *kinship_file_Output.txt* and *kinship_Output.txt*. The last keyword specifies the analysis parameter: *repetitions*. The text after the '=' are the keyword values.

### Keywords<a id="keywords-table"></a>
This is a list of OpenMendel keywords specific to Kinship. A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table). The names of keywords are *not* case sensitive. (The keyword values *may* be case sensitive.)


 Keyword          |   Default Value    | Allowed Values |  Short Description       
----------------      |  ----------------       |  ----------------      |  ----------------
kinship_output_file  |Kinship_Output_File.txt | User defined output file name | OpenMendel generated output file with table of kinship coefficients 
repetitions          |        1           |             integer            |  Repetitions for sharing statistics
xlinked_analysis  |  FALSE  |  TRUE, FALSE  |  Whether or not markers are on the X chromosome


### Data Files
Kinship requires a [Control file](https://openmendel.github.io/MendelBase.jl/#control-file), and a [Pedigree file](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Genotype data can be included in the Pedigree file, in which case a [Locus file](https://openmendel.github.io/MendelBase.jl/#locus-file) is required. Alternatively, genotype data can be provided in a [SNP data file](https://openmendel.github.io/MendelBase.jl/#snp-data-file), in which case a [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file) is required. OpenMendel will also accept [PLINK format](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml) FAM and BIM files. Details on the format and contents of the Control and data files can be found on the [MendelBase](https://openmendel.github.io/MendelBase.jl) documentation page. There are example data files in the Kinship [docs](https://github.com/OpenMendel/MendelKinship.jl/tree/master/docs) folder.

### Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelKinship

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in the control file Control_file.txt use the command:

     julia> Kinship("Control_file.txt")

*Note: The package is called* MendelKinship *but the analysis function is called simply* Kinship.

<!--- ### Interpreting the results
... --->

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
