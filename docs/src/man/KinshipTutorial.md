
# Tutorial of MendelKinship 


```julia
versioninfo() # system info for reproducibility
```

    Julia Version 1.5.0
    Commit 96786e22cc (2020-08-01 23:44 UTC)
    Platform Info:
      OS: macOS (x86_64-apple-darwin18.7.0)
      CPU: Intel(R) Core(TM) i9-9880H CPU @ 2.30GHz
      WORD_SIZE: 64
      LIBM: libopenlibm
      LLVM: libLLVM-9.0.1 (ORCJIT, skylake)


## About MendelKinship.jl and control files

This is a Julia implementation of the [Mendel software](https://academic.oup.com/bioinformatics/article/29/12/1568/292810) option 29, programmed in Fortran. `MendelKinship.jl` contains new features, but its usage largely resembles Mendel by employing *control files* - requiring users to make a file to specify options instead of interacting with Julia directly. This is not the most Julian way of programming or usage experience... but unfortuantely we lack manpower to correct it. 

## When to use MendelKinship

`MendelKinship.jl` is capable of calculating the theoretical kinship coefficient $\Phi_{ij}$ as long as a [valid pedigree structure](https://openmendel.github.io/MendelBase.jl/#pedigree-file) is provided. When SNP markers are available, `MendelKinship.jl` can also calculate empirical kinship coefficients using GRM, robust GRM, or Method of Moments methods (see [this paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20584) and [this paper](https://academic.oup.com/bioinformatics/article/26/22/2867/228512) for details). Here we recommend the Robust GRM or MoM (default) method because their estimates are more robust in the presence of rare alleles. 

`MendelKinship` can optionally compare the empirical kinship and theoretical kinships to check suspect pedigree structures and reveal hidden relatedness. It can also reveal sample mixed ups or other laboratory errors that can lead to inaccurate empirical kinships.  The result is saved in a table sorted in descending order. We optioanlly output 2 interactive plots that allow users to quickly pinpoint pairs with the greatest theoretical vs empirical deviance. 

## Examples data

Input data for this tutorial can be obtained on [our Github](https://github.com/OpenMendel/MendelKinship.jl/tree/master/data/documentation_data), which were originally derived from the 1000 genome project. They contain 85 people and 253141 SNPs, half of which have maf$< 0.05$. Using these founders' genotype, we simulated 127 extra people, resulting in 27 pedigrees and 212 people. Although the 85 individuals are treated as founders, they were actually somewhat related, and this is reflected in the kinship comparison in the 2nd example below. For more information on this dataset, please see Mendel's documentation example 29.4. 

## Using PLINK compressed file as input

MendelKinship additionally accepts [PLINK binary format](https://www.cog-genomics.org/plink2/formats#bed) as input, in which case the triplets (`data.bim`, `data.bed`, `data.fam`) must all be present. In this tutorial, there are no examples that uses these to import pedigree and SNP information. But if available, one can import the data by specifying the following in the control file:

`plink_input_basename = data` 

However, sometimes the .fam file contains non-unique person id (2nd column of .fam file) across different pedigrees, which is currently **not** permitted in MendelKinship. A person's id cannot be repeated in other pedigrees, even if it is contextually clear that they are different persons. This will be fixed in the near future.

## Analysis keywords available to users 

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

A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table)

## Example 1: Theoretical Kinship Coefficient Calculation 

### Step 1: Preparing the pedigree files:
Recall what is a [valid pedigree structure](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Note that we require a header line. The extension `.in` have no particular meaning. Let's examine (the first few lines of) such an example:


```julia
;head -10 "Ped29a.in"
```

    Pedigree,Person,Mother,Father,Sex,,,simTrait
      1       ,  16      ,          ,          ,  F       ,          ,  29.20564,
      1       ,  8228    ,          ,          ,  F       ,          ,  31.80179,
      1       ,  17008   ,          ,          ,  M       ,          ,  37.82143,
      1       ,  9218    ,  17008   ,  16      ,  M       ,          ,  35.08036,
      1       ,  3226    ,  9218    ,  8228    ,  F       ,          ,  28.32902,
      2       ,  29      ,          ,          ,  F       ,          ,  36.17929,
      2       ,  2294    ,          ,          ,  M       ,          ,  42.88099,
      2       ,  3416    ,          ,          ,  M       ,          ,  40.98316,
      2       ,  17893   ,  2294    ,  29      ,  F       ,          ,  35.55038,


### Step 2: Preparing the control file
A control file gives specific instructions to `MendelKinship`. To perform theoretical kinship calculation, an minimal control file looks like the following:


```julia
;cat "control_just_theoretical_29a.txt"
```

    #
    # Input and Output files.
    #
    pedigree_file = Ped29a.in
    #
    # Analysis parameters for Kinship option.
    #
    kinship_output_file = just_theoretical_output.txt

### Step 3: Run the analysis in Julia REPL or directly in notebook

We used the package Suppressor to hide warnings. They will be removed when we update `MendelKinship` to Julia version 1.0. However often informative warnings and/or MendelKinship messages will be printed, so it is best practice for new users to at least review the messages.


```julia
using MendelKinship
Kinship("control_just_theoretical_29a.txt")
```


<script>
// Immediately-invoked-function-expression to avoid global variables.
(function() {
    var warning_div = document.getElementById("webio-warning-11759167131154940107");
    var hide = function () {
        var script = document.getElementById("webio-setup-10996284096970835275");
        var parent = script && script.parentElement;
        var grandparent = parent && parent.parentElement;
        if (grandparent) {
            grandparent.style.display = "none";
        }
        warning_div.style.display = "none";
    };
    if (typeof Jupyter !== "undefined") {
        console.log("WebIO detected Jupyter notebook environment.");
        // Jupyter notebook.
        var extensions = (
            Jupyter
            && Jupyter.notebook.config.data
            && Jupyter.notebook.config.data.load_extensions
        );
        if (extensions && extensions["webio-jupyter-notebook"]) {
            // Extension already loaded.
            console.log("Jupyter WebIO nbextension detected; not loading ad-hoc.");
            hide();
            return;
        }
    } else if (window.location.pathname.includes("/lab")) {
        // Guessing JupyterLa
        console.log("Jupyter Lab detected; make sure the @webio/jupyter-lab-provider labextension is installed.");
        hide();
        return;
    }
})();

</script>
<p
    id="webio-warning-11759167131154940107"
    class="output_text output_stderr"
    style="padding: 1em; font-weight: bold;"
>
    Unable to load WebIO. Please make sure WebIO works for your Jupyter client.
    For troubleshooting, please see <a href="https://juliagizmos.github.io/WebIO.jl/latest/providers/ijulia/">
    the WebIO/IJulia documentation</a>.
    <!-- TODO: link to installation docs. -->
</p>



    env: node: No such file or directory


     
     
         Welcome to OpenMendel's
         Kinship analysis option
     
     
    Reading the data.
    
    The current working directory is "/Users/biona001/.julia/dev/MendelKinship/docs/src/man".
    
    Keywords modified by the user:
    
      control_file = control_just_theoretical_29a.txt
      kinship_output_file = just_theoretical_output.txt
      pedigree_file = Ped29a.in
     
     
    Analyzing the data.



<table class="data-frame"><thead><tr><th></th><th>Pedigree</th><th>Person1</th><th>Person2</th><th>Kinship</th><th>Delta7</th><th>delta1</th><th>delta2</th><th>delta3</th><th>delta4</th></tr><tr><th></th><th>String</th><th>String</th><th>String</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th><th>Float64</th></tr></thead><tbody><p>1,833 rows × 14 columns (omitted printing of 5 columns)</p><tr><th>1</th><td>1</td><td>16</td><td>16</td><td>0.5</td><td>1.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>2</th><td>1</td><td>16</td><td>17008</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>3</th><td>1</td><td>16</td><td>3226</td><td>0.125</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>4</th><td>1</td><td>16</td><td>8228</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>5</th><td>1</td><td>16</td><td>9218</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>6</th><td>1</td><td>17008</td><td>17008</td><td>0.5</td><td>1.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>7</th><td>1</td><td>17008</td><td>3226</td><td>0.125</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>8</th><td>1</td><td>17008</td><td>9218</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>9</th><td>1</td><td>3226</td><td>3226</td><td>0.5</td><td>1.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>10</th><td>1</td><td>8228</td><td>17008</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>11</th><td>1</td><td>8228</td><td>3226</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>12</th><td>1</td><td>8228</td><td>8228</td><td>0.5</td><td>1.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>13</th><td>1</td><td>8228</td><td>9218</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>14</th><td>1</td><td>9218</td><td>3226</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>15</th><td>1</td><td>9218</td><td>9218</td><td>0.5</td><td>1.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>16</th><td>2</td><td>14695</td><td>14695</td><td>0.5</td><td>1.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>17</th><td>2</td><td>14695</td><td>17893</td><td>0.25</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>18</th><td>2</td><td>14695</td><td>6952</td><td>0.125</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>19</th><td>2</td><td>17893</td><td>17893</td><td>0.5</td><td>1.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>20</th><td>2</td><td>17893</td><td>6952</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>21</th><td>2</td><td>2294</td><td>14695</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>22</th><td>2</td><td>2294</td><td>17893</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>23</th><td>2</td><td>2294</td><td>2294</td><td>0.5</td><td>1.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>24</th><td>2</td><td>2294</td><td>3416</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>25</th><td>2</td><td>2294</td><td>3916</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>26</th><td>2</td><td>2294</td><td>6790</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>27</th><td>2</td><td>2294</td><td>6952</td><td>0.125</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>28</th><td>2</td><td>29</td><td>14695</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>29</th><td>2</td><td>29</td><td>17893</td><td>0.25</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>30</th><td>2</td><td>29</td><td>2294</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table>


    
     
     
    Mendel's analysis is finished.
    


### Step 4: Interpreting the result

`MendelKinship` should have generated the file`just_theoretical_output.txt` in your local directory. One can directly open the file, or import into the Julia environment for ease of manipulation using the DataFrames package. The fourth column contains the desired theoretical kinship coefficient. The 5th column contains the (deterministically) estimated Delta7 matrix. The 6th through the 14 columns contain the (stochastically) estimated Jacquard's 9 identity coefficients.

## Example 2: Compare theoretical/empirical kinship values

When both pedigree structure and *complete* SNP information are available, we can compare theoretical/empirical kinship coefficients. In practice, however, we often have individuals without genotype information, but nevertheless must be included in the pedigree structure. `MendelKinship` does not handle this situation yet, but an analysis option that supports these data is being developed. For now you can impute genotypes but keep in mind that the relationship comparison for these individuals who lack all genotype information will not be meaningful.   

### Step 1: Prepare pedigree file and SNP data file

The pedigree file is the same as the pedigree file in the previous example. The SNP definition file requires a header row, and should have approprietely placed commas. It may be informative to compare the following SNP definition file with the original "SNP_def29a.in" in Mendel Option 29a. 


```julia
;head -10 "SNP_def29a_converted.txt"
```

    Locus,Chromosome,Basepairs,Allele1,Allele2
    rs3020701,19,90974,1,2
    rs56343121,19,91106,1,2
    rs143501051,19,93542,1,2
    rs56182540,19,95981,1,2
    rs7260412,19,105021,1,2
    rs11669393,19,107866,1,2
    rs181646587,19,107894,1,2
    rs8106297,19,107958,1,2
    rs8106302,19,107962,1,2


#### Non binary PLINK users

The SNP data files in this case must be stored in PLINK BED file in SNP-major format, with an accompanying SNP definition file. For an explanation of what these are, see [MendelBase documentation](https://openmendel.github.io/MendelBase.jl/).

#### Binary PLINK file users

If your have "data.bim", "data.bed", "data.fam" (i.e. the 3 triplet of PLINK files), then you can replace the 3 fields `snpdata_file`, `snpdefinition_file`, and `pedigree_file` in the next step with just 1 field:

`plink_input_basename = data`.

### Step 2: Preparing control file

The following control file tells MendelKinship to compare theoretical kinship and empirical kinship, and output 2 interactive plots stored in .html format. 


```julia
;cat "control_compare_29a.txt"
```

    #
    # Input and Output files.
    #
    snpdata_file = SNP_data29a.bed
    snpdefinition_file = SNP_def29a_converted.txt
    pedigree_file = Ped29a.in
    #
    # Analysis parameters for Kinship option.
    #
    compare_kinships = true
    kinship_plot = kinship_plot
    z_score_plot = z_score_plot

### Step 3: Running the analysis


```julia
using MendelKinship
Kinship("control_compare_29a.txt")
```

     
     
         Welcome to OpenMendel's
         Kinship analysis option
     
     
    Reading the data.
    
    The current working directory is "/Users/biona001/.julia/dev/MendelKinship/docs/src/man".
    
    Keywords modified by the user:
    
      compare_kinships = true
      control_file = control_compare_29a.txt
      kinship_plot = kinship_plot
      pedigree_file = Ped29a.in
      snpdata_file = SNP_data29a.bed
      snpdefinition_file = SNP_def29a_converted.txt
      z_score_plot = z_score_plot
     
     
    Analyzing the data.
    
    Kinship plot saved.
    Fisher's plot saved.


<table class="data-frame"><thead><tr><th></th><th>Pedigree1</th><th>Pedigree2</th><th>Person1</th><th>Person2</th><th>theoretical_kinship</th><th>empiric_kinship</th><th>fishers_zscore</th></tr><tr><th></th><th>String</th><th>String</th><th>String</th><th>String</th><th>Float64</th><th>Float64</th><th>Float64</th></tr></thead><tbody><p>22,578 rows × 7 columns</p><tr><th>1</th><td>14</td><td>14</td><td>26732</td><td>264</td><td>0.0</td><td>0.109552</td><td>5.31942</td></tr><tr><th>2</th><td>31</td><td>31</td><td>15884</td><td>19770</td><td>0.25</td><td>0.150364</td><td>-4.2355</td></tr><tr><th>3</th><td>23</td><td>23</td><td>9943</td><td>392</td><td>0.125</td><td>0.0225133</td><td>-4.20155</td></tr><tr><th>4</th><td>25</td><td>14</td><td>22041</td><td>16636</td><td>0.0</td><td>0.0969715</td><td>4.75137</td></tr><tr><th>5</th><td>25</td><td>25</td><td>11822</td><td>24192</td><td>0.25</td><td>0.159229</td><td>-3.82975</td></tr><tr><th>6</th><td>14</td><td>14</td><td>25732</td><td>264</td><td>0.125</td><td>0.216622</td><td>4.62515</td></tr><tr><th>7</th><td>25</td><td>25</td><td>3012</td><td>3016</td><td>0.125</td><td>0.213888</td><td>4.4971</td></tr><tr><th>8</th><td>25</td><td>17</td><td>23404</td><td>12004</td><td>0.0</td><td>-0.0896437</td><td>-3.60943</td></tr><tr><th>9</th><td>25</td><td>23</td><td>23404</td><td>19279</td><td>0.0</td><td>-0.0877967</td><td>-3.52627</td></tr><tr><th>10</th><td>17</td><td>14</td><td>26857</td><td>264</td><td>0.0</td><td>0.0859953</td><td>4.25691</td></tr><tr><th>11</th><td>10040</td><td>8</td><td>234</td><td>5226</td><td>0.0</td><td>-0.0859709</td><td>-3.44408</td></tr><tr><th>12</th><td>23</td><td>19</td><td>743</td><td>8344</td><td>0.0</td><td>-0.0856543</td><td>-3.42983</td></tr><tr><th>13</th><td>19</td><td>19</td><td>22375</td><td>720</td><td>0.125</td><td>0.0405078</td><td>-3.39689</td></tr><tr><th>14</th><td>25</td><td>23</td><td>23404</td><td>743</td><td>0.0</td><td>-0.0831952</td><td>-3.3192</td></tr><tr><th>15</th><td>10040</td><td>8</td><td>234</td><td>13234</td><td>0.0</td><td>-0.0805884</td><td>-3.20196</td></tr><tr><th>16</th><td>23</td><td>8</td><td>3121</td><td>19099</td><td>0.0</td><td>-0.0792375</td><td>-3.14122</td></tr><tr><th>17</th><td>10040</td><td>8</td><td>234</td><td>7395</td><td>0.0</td><td>-0.079153</td><td>-3.13743</td></tr><tr><th>18</th><td>10040</td><td>8</td><td>234</td><td>19099</td><td>0.0</td><td>-0.0789209</td><td>-3.12699</td></tr><tr><th>19</th><td>17</td><td>17</td><td>908</td><td>10418</td><td>0.25</td><td>0.174596</td><td>-3.12361</td></tr><tr><th>20</th><td>31</td><td>23</td><td>19770</td><td>3760</td><td>0.0</td><td>-0.0783826</td><td>-3.1028</td></tr><tr><th>21</th><td>31</td><td>31</td><td>152</td><td>19770</td><td>0.125</td><td>0.0472729</td><td>-3.0941</td></tr><tr><th>22</th><td>14</td><td>14</td><td>25732</td><td>25732</td><td>0.5</td><td>0.55628</td><td>3.89575</td></tr><tr><th>23</th><td>10040</td><td>11</td><td>234</td><td>7670</td><td>0.0</td><td>-0.0772533</td><td>-3.05204</td></tr><tr><th>24</th><td>25</td><td>17</td><td>24192</td><td>12004</td><td>0.0</td><td>-0.0769367</td><td>-3.03781</td></tr><tr><th>25</th><td>25</td><td>17</td><td>23404</td><td>19113</td><td>0.0</td><td>-0.0769156</td><td>-3.03687</td></tr><tr><th>26</th><td>23</td><td>19</td><td>743</td><td>23699</td><td>0.0</td><td>-0.0768523</td><td>-3.03402</td></tr><tr><th>27</th><td>17</td><td>8</td><td>12004</td><td>5226</td><td>0.0</td><td>-0.0767573</td><td>-3.02975</td></tr><tr><th>28</th><td>23</td><td>19</td><td>3760</td><td>8344</td><td>0.0</td><td>-0.0763879</td><td>-3.01315</td></tr><tr><th>29</th><td>10029</td><td>17</td><td>434</td><td>12004</td><td>0.0</td><td>-0.0760713</td><td>-2.99893</td></tr><tr><th>30</th><td>10040</td><td>25</td><td>234</td><td>3012</td><td>0.0</td><td>-0.0760713</td><td>-2.99893</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table>


    
     
     
    Mendel's analysis is finished.
    


### Step 4: Interpreting the Result

Founders which have 0 theoretical kinships often exhibit a non-zero empirical kinship. In the first row, person 26732 and 264 have 0 theoretical kinship but their empirical kinship is pretty close to 0.125 = 1/8. That is, these 2 people which we initially thought are unrelated, may be half siblings, grandparent-grandchild, or an avuncular pair. On the otherhand, the 8th row has a founder pair that has a $-0.08$ kinship (i.e. they are very *un*related), suggesting that the standard deviation of the moments estimator may have a wide spread. There may also have been a sample mix up. Another explanation is that we are only using one chromosome's worth of data and so the estimates of kinship may be imprecise. 

## Interactive Plots and Tables

`MendelKinship` automatically generates 2 figures and 1 table to allow the user to easily compare theoretical and empirical kinship, detect outliers, and observe skewnesses in distribution. Figures are saved in `.html` format to enable interactive sessions. To summarize, 

+ The table containing all the pairwise kinship and theoretical comparisons is stored in `kinship_file_output.txt`. The table is sorted in descending order of the largest deviance between the theoretical and empiric kinship. The last column lists the [Fisher's Z statistic](https://en.wikipedia.org/wiki/Fisher_transformation) (i.e. the number of standard deviations away from mean). 
    
+ The 2 plots are stored in .html format, which should be automatically be generated in your directory. These figures can be examined interactively via jupyter notebook, as demonstrated below, or opened directly via the browser.

### Generated Interactive Plots part 1:

The first interactive plot allows user to quickly identify which pairs of persons have an empirical kinship most deviated from their expected (theoretical) kinship. The midpoint is placed as an orange dot for interpretability. As an example, the first row in the table above is the highest point on the left most spread. Careful readers might observe that there is a wider spread on those with 0 expected theoretical kinship. This is expected, because most people are not related to each other, so we are making many more comparisons that have 0 expected kinship. 


```julia
using MendelKinship, PlotlyJS, CSV

#import the files created from the previous example
result = CSV.read("kinship_file_output.txt")
name = Vector{String}(undef, size(result, 1))

# label the data points according to the persons names
for i in 1:length(name)
    name[i] = "Person1=" * string(result[i, 3]) * ", " * "Person2=" * string(result[i, 4])
end

#create interactive graph
function compare_kinship_plot()
    trace1 = scatter(;x=result[:theoretical_kinship], 
        y=result[:empiric_kinship], mode="markers", 
        name="empiric kinship", text=name)
    
    trace2 = scatter(;x=[1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 0.0],
        y=[1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 0.0], 
        mode="markers", name="marker for midpoint")
        
    layout = Layout(;title="Compare empiric vs theoretical kinship",hovermode="closest", 
        xaxis=attr(title="Theoretical kinship (θ)", showgrid=false, zeroline=false),
        yaxis=attr(title="Empiric Kinship", zeroline=false))
    
    data = [trace1, trace2]
    plot(data, layout)
end
compare_kinship_plot()
```




<div
    class="webio-mountpoint"
    data-webio-mountpoint="13921852083029915973"
>
    <script>
    if (window.require && require.defined && require.defined("nbextensions/webio-jupyter-notebook")) {
        console.log("Jupyter WebIO extension detected, not mounting.");
    } else if (window.WebIO) {
        WebIO.mount(
            document.querySelector('[data-webio-mountpoint="13921852083029915973"]'),
            window,
        );
    } else {
        document
            .querySelector('[data-webio-mountpoint="13921852083029915973"]')
            .innerHTML = (
                '<div style="padding: 1em; background-color: #f8d6da; border: 1px solid #f5c6cb">' +
                '<p><strong>WebIO not detected.</strong></p>' +
                '<p>Please read ' +
                '<a href="https://juliagizmos.github.io/WebIO.jl/latest/troubleshooting/not-detected/" target="_blank">the troubleshooting guide</a> ' +
                'for more information on how to resolve this issue.</p>' +
                '<p><a href="https://juliagizmos.github.io/WebIO.jl/latest/troubleshooting/not-detected/" target="_blank">https://juliagizmos.github.io/WebIO.jl/latest/troubleshooting/not-detected/</a></p>' +
                '</div>'
            );
    }
    </script>
</div>




### Generated Interactive Plots part 2:

After comparing the theoretical and empirical kinships, as in the previous graph or through the outputted table directly, often one may wonder whether the observed differences between the two statistics are significantly different. As explained in our main OpenMendel paper (section 7), this difference can be precisely quantified by the Fisher's z transformation, which should give us samples from a standard normal distribution $N(0, 1)$. We ploted this statistic in plot 2, and at first glance, the distribution is approximately normal. In Julia, we can easily verify this by computing some summary statistics:


```julia
function fishers_transform()
    trace1 = histogram(x=result[:fishers_zscore], text=name)
    data = [trace1]
    
    layout = Layout(barmode="overlay", 
        title="Z-score plot for Fisher's statistic",
        xaxis=attr(title="Standard deviations"),
        yaxis=attr(title="count"))
    
    plot(data, layout)
end
fishers_transform()
```




<div
    class="webio-mountpoint"
    data-webio-mountpoint="3685400258333464949"
>
    <script>
    if (window.require && require.defined && require.defined("nbextensions/webio-jupyter-notebook")) {
        console.log("Jupyter WebIO extension detected, not mounting.");
    } else if (window.WebIO) {
        WebIO.mount(
            document.querySelector('[data-webio-mountpoint="3685400258333464949"]'),
            window,
        );
    } else {
        document
            .querySelector('[data-webio-mountpoint="3685400258333464949"]')
            .innerHTML = (
                '<div style="padding: 1em; background-color: #f8d6da; border: 1px solid #f5c6cb">' +
                '<p><strong>WebIO not detected.</strong></p>' +
                '<p>Please read ' +
                '<a href="https://juliagizmos.github.io/WebIO.jl/latest/troubleshooting/not-detected/" target="_blank">the troubleshooting guide</a> ' +
                'for more information on how to resolve this issue.</p>' +
                '<p><a href="https://juliagizmos.github.io/WebIO.jl/latest/troubleshooting/not-detected/" target="_blank">https://juliagizmos.github.io/WebIO.jl/latest/troubleshooting/not-detected/</a></p>' +
                '</div>'
            );
    }
    </script>
</div>




### Compute mean and variance

We can verify that the Fisher's statistic is approximately normal by checking its 1~4th moments:


```julia
using Statistics, StatsBase

my_zscore = convert(Vector{Float64}, result[:fishers_zscore])
mean(my_zscore), var(my_zscore)
```




    (-1.2084702388691578e-16, 1.0)



### Compute skewness and excess kurtosis


```julia
skewness(my_zscore), kurtosis(my_zscore)
```




    (0.0923976421835524, 0.10657222736224581)



## Conclusions
MendelKinship provides a rapid way to calculate the theoretical kinship, which requires accurate pedigrees. Calculation of empirical kinships is also possible if genotypes at multiple markers are availble. Further, it can compare these kinships using Fisher's Z statistic when both the pedigrees and markers are available. 