![GitHub Logo](logo.png)

#### <v1.0.1>

## About

**Manticore** addresses the subgenomic intermixing state of a hybrid, allopolyploid genome. The user provides a reference FASTA file for the hybrid species and short illumina reads from its two parental species. Manticore studies the coverage profiles obtained from each parent and extracts a series of metrics in intervals of adjustable size.

The provided metrics indicate how conserved the subgenomes are when compared to their parental sequences. They also indicate how reshuffled the two subgenomes are, a metric we define as **subgenomic intermixing** and we address extensively in the associated publication (**submitted**).

The user may use **Manticore** as collateral information when studying the genome of a hybrid species. In fact, the only necessary data for the tool to run are a hybrid genome in FASTA format and parental reads in FASTQ format. Especially when working with complex genomes such as those of plants, we recommend the usage of our tool to address the subgenomic conservation and intermixing state.

**Manticore** produces all its output inside the user-declared `--output-dir`. Inside this folder, six subdirectories are found, numerically ordered by time of creation. Folders `1_` to `4_` contain only intermediate files, and can be disregarded unless the user wants to address anything specifically related to **bam** or **depth** files. The results are provided in the `5_` and `6_` subdirectories. The tables and metrics are provided within the `5_tables` subdirectory. The plots are found within the `6_plots` subdirectory. The plots are provided both in *svg* and *png* format, together with the *R* scripts used to generate them. This allows the users to tweak the plot to their liking.

## Requirements

Manticore depends on a series of external resources, which are summarized in the following table:

| Resource | Type | Version | Source |
|-|-|-|-|
| Python | Language interpreter | ≥ 3.5 | https://www.python.org/downloads |
| R | Language interpreter | ≥ 3.6.0 | https://www.r-project.org |
| HISAT2 | Mapping program | ≥ 2.1.0 | https://github.com/DaehwanKimLab/hisat2 |
| Biopython | Python package | ≥ 1.73 | python3.6 -m pip install biopython |
| pysam | Python package | ≥ 0.15.3 | python3.6 -m pip install pysam |
| numpy | Python package | ≥ 1.17.3 | python3.6 -m pip install numpy |
| pandas | Python package | ≥ 1.0.3 | python3.6 -m pip install pandas |
| dask | Python package | ≥ 2.12.0 | python3.6 -m pip install "dask[complete]" |
| ggplot2 | R package | ≥ 3.3.2 | https://ggplot2.tidyverse.org/ |

</bg>

*HISAT2* (Kim et al., 2015) is used as an external software for the read mapping, and *R* is used for the plots generation. These programs are automatically searched in the `$PATH`, and if not found, an error is raised. A user can however specify the path to the executable of each one of them (see command line options), which will override the `$PATH` search.

Manticore is based on *python3*. Besides the python modules listed in the table above, Manticore also depends on the following built-in python modules, which shouldn't require an installation as they come together with a *python3* installation: `argparse, itertools, math, multiprocessing, operator, os, os.path, shutil, subprocess, sys, time`. Manticore was developed on a high-performance computing (HPC) cluster. Depending on the read data set used, its RAM usage might be quite high. While it is absolutely possible to run it on a common laptop, with complex or large genomes we recommend the average user to have at least 30 Gb RAM available. The core analysis of Manticore is performed making use of thread parallelization. Hence, the more cpu cores will be available, the quicker the results will come.

## Installation

No installation step is required from the user, apart from `git clone git@github.com:MatteoSchiavinato/manticore.git`. The user shall keep in mind that the Manticore main script (`./manticore`) uses the `./src/analyze-windows.py` script. Hence, the `./src` subdirectory must always be where the main script is.

## Usage

The program has many options, but only few of them are mandatory. If the user runs the program with minimal parameters, it will just require:

1) a reference FASTA sequence
2) FASTQ paired-end reads from the parents
3) an output directory

The whole analysis will be run within such folder, making its management easy from the user at the end of the analysis and keeping the results tidy. We do recommend to set `--max-mem` and `--threads` for a shorter computation time of the analysis.

### Basic run

By running `./manticore -h` the user accesses a quick help that summarizes the options available. This help is intended for the user that already is comfortable with the available options. If this is not the case, an extended help can be printed by running `./manticore --help`.

The quick help shows the mandatory options first, which are the only ones actually required to run the program (i.e. all the other options are fine-tuning). An example command is provided here, with explanation below:

```
Usage:
manticore \
--species-name STRING \
--reads-type PE \
--reads PATH1,PATH2 PATH3,PATH4 \
        (parentA)   (parentB)
--names STRING STRING \
        (A)    (B)
--reference PATH \
--output-dir PATH \
[ ... other options ... ]
```

It is important to notice the structure of the`--reads` argument, as it is crucial for a successful run. This argument is a SPACE-separated list. Each element in the list has read file 1 and 2 from a parent of the hybrid, separated by a comma (`parent_A.read_1.fq,parent_A.read_2.fq`). In this example, parent A and parent B are declared in this order (i.e. A and B). The `--names` argument assigns a “name” to each parent using the same order. This means that “parent A” will be referred to with the first STRING declared in `--names`, and the second will be used for parent B. The names are chosen by the user and, as long as they don't contain whitespaces or special characters (e.g. `$^/...`), they won't be an issue.

### Resuming a stopped run

If your run happens to crash at a certain point, Manticore has a built-in system to resume it. A series of files with the `*.done` extension are created during the analysis: these files help the program understand which steps have already been performed. When running the program with `--output-dir`, the program checks this folder for existing `*.done` files and skips all the steps for which a `*.done` file is already present.

The removal of a `*.done` file causes the program to re-run the associated step; a user can also specify the `--restart` option to override all the `*.done` files and start from scratch. The `*.done` files shall not be moved or renamed, as in that case they won't be found anymore.

The more experienced users may *trick* the program by removing the `*.done` file of the step that they want to reproduce. For example, the program performed the mapping and the filtering but then crashed at the sorting step. The users want to resume that run, but decided meanwhile that they want to repeat the filtering too, changing filtering parameters with the `--samtools-filters` option. Hence, they delete the associated `*.done` file within the mapping subdirectory and this will lure Manticore into repeating the filtering step.

### Results - tables

Manticore produces two types of results. The sub-folder `5_tables` will contain several TAB-separated tables with the computed information. The `combined.results.txt` file contains all the results produced for each window. Its content is structured like this:

```
Scaffold       W_start   W_end     Feature  Real_length  Jaccard  Subgenome  Cov_pos  Frac_pos  Mean_cov  Union   Intersection  Uncovered
chrA01         1         500000    CDS      97237        0.0      Genome_A   92020    95.28     23.31     500000  0             0
chrA01         1         500000    CDS      97237        0.0      Genome_C   14961    14.41     11.94     500000  0             0
chrA01         1         500000    nonrep   412677       0.1      Genome_A   369000   89.08     23.98     500000  50000         0
chrA01         1         500000    nonrep   412677       0.1      Genome_C   53194    13.35     21.51     500000  50000         0
chrA01         500001    1000000   CDS      139951       0.0      Genome_A   133987   95.46     23.03     500000  0             0
chrA01         500001    1000000   CDS      139951       0.0      Genome_C   14591    10.76     10.75     500000  0             0
chrA01         500001    1000000   nonrep   449655       0.0      Genome_A   408979   90.92     22.06     500000  0             0
chrA01         500001    1000000   nonrep   449655       0.0      Genome_C   45169    10.17     12.53     500000  0             0
chrA01         1000001   1500000   CDS      143992       0.0      Genome_A   137496   95.34     21.76     500000  0             0
```

Each line represents a window in a certain sequence (columns `Scaffold`, `W_start`, `W_end`). Multiple lines can point at the same window if more features (i.e. BED files) are passed (column `Feature`).

The `Real_length` column indicates the number of positions within said window that were overlapping the BED file indicated in column `Feature`. If no BED file was passed, real length shall be the same length as `W_end - W_start + 1`.

The `Jaccard` column contains the computed Jaccard index between the coverage profiles of the two parents within the window; details on this computation can be found in the associated publication (**submitted**). The `Subgenome` column shows which parental coverage profile does the line correspond to.

The `Cov_pos`, `Frac_pos` and `Mean_cov` columns show, respectively: the number of covered positions; the fraction of covered positions (in %), the mean coverage of the window.

The `Union`, `Intersection` and `Uncovered` columns show, respectively: the union of the positions covered by both parents; the intersection of those positions; the positions left uncovered. These last columns contain values that are multiples of `--window-size` divided by `--n-breaks`. This length is called **break length**.

The `5_tables` folder also contains other files. The `combined.results.txt` output file is also present in a split version, in as many sub-files as there are Subgenomes and Features. If the passed subgenomes are `(Genome_A, Genome_C)`, and the passed features are `(CDS, nonrep)`, there will be 2 x 2 = 4 files named as `<species_name>.<subgenome>.<feature>.txt`. These files contain the same information as `combined.results.txt`, but partitioned by Subgenome and Feature. The `<feature>.relative_lengths.table` contains the information in base pairs and in percentage about the genome and the amount of it that is represented in the windows. The `<species_name>.coverage_metrics.txt` file shows other information by Feature (1st column) and Subgenome (2nd column). When a dash (-) is present in the 2nd column, the metric refers to the entire genome and not to one specific subgenome. Here is a quick overview of the metrics:

| Metric | Meaning |
| - | - |
| Total_length | Total base pairs contained in the analysed windows |
| Featured_length | Total base pairs annotated as the feature in Column 1 |
| Union | Union of the positions covered by Parent 1 and 2 (logical OR) |
| Intersection | Intersection of the positions covered by Parent 1 and 2 (logical AND) |
| Uncovered | Positions where neither parent produced coverage |

The file also contains some subgenome-specific metrics. The subgenome is indicated in column 2. Here is a quick overview of the metrics:

| Metric | Meaning |
| - | - |
| Total | Sum of the positions covered by this parent |
| Assigned\_windows | Number of windows that could be assigned to this parent (based on `--max-jacc-uniq`) |
| Assigned\_window\_pos | Base pairs assigned to this parent (based on "Assigned\_windows" * `--window_size`) |
| Unique | Base pairs uniquely covered by this parent |
| Mean_coverage | Mean coverage of the windows assigned to this parent |
| Mean\_cov_frac | Mean covered fraction (%) of the windows assigned to this parent |
| Mean\_cov\_frac_0-50 | ... of the windows covered <\ 50% |
| Mean\_cov\_frac_50-100 | ... of the windows covered from 50% to 100% |

### Results - plots

Manticore also produces plots for the direct interpretation of the results from the user. \[to be continued\]

### Command-line options

The following is a list of all the command line options of Manticore:

**[INPUT OPTIONS]**

`--species-name`<br/>
This word will be used as a prefix for many files throughout the analysis. It is suggested to use a single word that corresponds to the hybrid species that is being studied.
[mandatory]

`--reads`<br/>
SPACE-separated list of paired-end read files. Each SPACE should separate reads from different parents, while the two read files (R1 and R2) of the same parent should be separated by a COMMA (without a following space). An example is provided here: Par_1.R1.fq,Par_1.R2.fq Par_2.R1.fq,Par_2.R2.fq
In case of single-end reads, ignore the comma separation and only provide the two parental file separated by a space.
[mandatory]

`--reads-type`<br/>
Specify either 'SE' or 'PE' depending on if your reads are single-end or paired-end, respectively.
[PE]

`--names`<br/>
SPACE-separated list of names to be attributed to the parental reads listed in `--reads`. The order of the names has to be the same of `--reads` (i.e. first name will be assigned to first pair of reads, and so on).
[mandatory]

`--reference`<br/>
FASTA file where to assess genomic intermixing. This is the FASTA file corresponding to the genome assembly of the hybrid species that is being studied. It doesn't have to be a chromosome-level assembly, although fragmented assemblies slow down the analysis sensibly.
[mandatory]

`--output-dir`<br/>
All files will be produced within this directory. It is suggested not to use the current (“.”), as many folders will be produced and they may unluckily overlap with folders that are already existing in current (“.”).
[mandatory]

`--region-beds`<br/>  
SPACE-separated list of BED files with regions that have to be considered. For each BED file (and the corresponding name specified in `--region-names`) a separate line of analysis will be run, eventually leading to a separate set of plots and tables. If no BED file is specified, the program will consider all the genome.
[-]

`--region-names`<br/>
SPACE-separated list of names to associate to BED files. We suggest to use names that define the type of data they contain, as these names will be used in the plots. Each name should be composed of one single word.
[-]


**[COVERAGE ANALYSIS]**

`--window-size`<br/>
Size of the windows on which to study intermixing. Larger windows will inevitably lead to more intermixing observed, and therefore to larger jaccard indexes. The opposite can be said for shorter windows.
[500000]

`--n-breaks`<br/>
Number of window breaks to compute Jaccard index from. `--window-size` divided by `--n-breaks` has to return an integer (Hint: the smaller, the more likely intermixing is observed)
[10]

`--min-feat-length`<br/>
Minimum length of annotated features within windows of size `--window-size`. This applies separately to all annotations provided with `--region-beds`. Each window must include this many positions annotated in BED files (if any); if not, the window will be assigned an 'NA' jaccard index. If no BED files are provided, this option just needs to be smaller than `--window-size`.
[1000]

`--min-frac-pos`<br/>
Minimum fraction of positions that have to be covered within the `--min-feat-length` of a window. This applies separately to all annotations provided with `--region-beds`. The union of the positions covered by all parental reads must cover at least this fraction of the positions defined by BED files (if any) in any window. If not, the window will be assigned an 'NA' jaccard index. This value is a float from 0 to 1.
[0.1]

`--min-cov-pos`<br/>
Minimum number of positions that have to be covered within the `--min-feat-length` of a window. This applies separately to all annotations provided with `--region-beds`. The union of the positions covered by all parental reads must sum up to at least this value in any window. If not, the window will be assigned an 'NA' jaccard index. This value is an integer.
[1000]

`--min-cov`<br/>
Minimum position coverage to consider a position in the analysis. Raising this value eventually shifts the jaccard index distribution towards 0, as there will be less overlap between the parental read coverage profiles. At the same time, raising this value reduces false positives.
[1]


**[TABLES & PLOTTING]**

`--max-jacc-uniq`<br/>
At the end of the analysis, Manticore goes through the results obtained for each window to assign each window to a parental progenitor. The assignment is performed only if the Jaccard index is lower or equal to the value specified with this option. This value is a float ranging from 0 to 1. Lower values will return less assigned positions (i.e. lower values are more stringent).
[0.5]

`--max-plot-cov`<br/>
The final plots will be limited to this maximum coverage. This value is only used at the stage of plotting.
[100]


**[MISCELLANEOUS]**

`--max-mem`<br/>
Maximum memory that can be used by the program (use only '#G'). N threads will get '--max-mem / N' memory each. The program will not use more memory than this. Note that specifying the maximum available memory of your hardware might generate pysam-related issues. A larger `--max-mem` will likely speed up the analysis. In HPC clusters, we suggest to not use more than 50G, as they won't likely be needed.
[4G]

`--threads`<br/>
Number of parallel threads.
[4]

`--cleanup`<br/>
Heavy intermediate files are deleted when the program has finished running. This means that the raw (i.e. unfiltered) mapping files, as well as the unfiltered depth files will be removed from their corresponding directories. This option will help the users with the need to control their memory usage. However, since the files are removed only at the end, one must be aware that this option does not reduce the memory usage during the run, rather just afterwards.
[off]

`--filter-reference`<br/>
The reference file is filtered, keeping only sequences longer than the specified `--window-size`. This option has both pros and cons: parental reads originating from regions that are orthologous to the excluded sequences won't find their mapping target; however, there will be less sequences to analyse and this will speed up the analysis.
[off]

`--isize-read-num`<br/>
Number of read pairs from which to estimate the insert size distribution. This number of reads will be mapped on the reference and used to draw an insert size distribution plot, which will then be used to define the peak insert size and the range of accepted insert sizes through `--isize-dist-width`.
This option is ignored in case of single-end reads.
[10000]

`--isize-dist-width`<br/>
Width of the allowed insert size range when mapping. The `--isize-read-num` option returns a peak insert size P. If a value of **N** is specified here, then this value is used in the following way: the -I and -X parameters of *HISAT2*, which control the insert size range, will be `P-N` and `P+N`, respectively. A wider width allows for more read pairs to be accepted as valid, but at the same time increases the chance of accepting wrongly mapped reads.
This option is ignored in case of single-end reads.
[150]

`--hisat2-path`<br/>
Path to the *HISAT2* executable (only specify if not present in the `$PATH`).
[hisat2]

`--hisat2-map-pars`<br/>
Mapping parameters to be passed to *HISAT2*. For a detailed description of these parameters, consult the HISAT2 manual directly ( https://ccb.jhu.edu/software/hisat2/manual.shtml ) directly.
[-k 5 --score-min L,0.0,-0.6 --mp 6,2 --rdg 5,3 --rfg 5,3 --no-softclip --no-spliced-alignment]

`--samtools-filters`<br/>
Arguments to pass to `samtools view`. Manticore is designed to only handle flag inclusion or exclusion in this field, that is, with `-F` and `-f` only. Using the other samtools arguments in here is *not* recommended. The default parameters exclude secondary alignments (`-F 0x0100`) and records of unmapped reads (`-F 0x4`). An optional parameter that a user might want to include is to retain only proper pairs (`-f 0x2`). We stress that this option sensibly reduces the amount of retained read pairs.
[-F 0x0100 -F 0x4]

`--rscript-path`<br/>
Path to the `Rscript` executable (only specify if not present in the `$PATH`).
[Rscript]

`--version`<br/>
Print program name, version and exit.
[-]

`--restart`<br/>
Ignore any existing `*.done` file and restart the analysis (each directory contains a `*.done` file that signals to skip the step).
[off]
