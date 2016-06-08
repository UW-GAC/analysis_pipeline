# TOPMed analysis pipeline

## Setup

1. Install R packages and dependencies from Bioconductor
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("SeqVarTools", "SNPRelate", "GENESIS", "argparser"))
```
2. Install updated GENESIS from github
```{r}
library(devtools)
install_github("smgogarten/GENESIS")
```
3. Install TopmedPipeline R package
```
R CMD INSTALL TopmedPipeline
```

## Basic outline

Each script in the `R` directory takes a config file with parameters. Look at the beginning of each script for parameter lists. Some parameters are required; others are optional with default values.
Example config files are in `testdata`.

Some scripts can be run in parallel by chromosome. For these scripts, the chromosome number is the second argument. If running in parallel, include a space in file names in the config file where chromosome should be inserted, e.g.,
```
gds_file "1KG_phase1_release_v3_chr .gds"
```

Nearly all scripts require a GDS file in SeqArray format. Phenotype files should be an AnnotatedDataFrame saved in an RData file. See `?AnnotatedDataFrame` or the SeqVarTools documentation for details. Example files are provided in `testdata`.

Python scripts are provided to run multi-step analyses on a compute
cluster or cloud environment. The `submitJob` function in `TopmedPipeline.py` is written
for a Sun Grid Engine (SGE) cluster, but may be modified for other
environments. These scripts require a config argument `out_prefix` in
addition to the arguments for each R script called. Some input and output
file name parameters are overridden by the scripts in order to link
jobs together.


## Conversion to GDS

`vcf2gds.py`

1. `vcf2gds.R`
2. `merge_gds.R`
3. `unique_variant_ids.R`

Step 1 converts VCF files (one per chromosome) into GDS files,
discarding non-genotype FORMAT fields. Step 2 combines these files
into a single GDS file, which is needed for whole-genome analyses such
as relatedness and population structure. The single-chromosome files
are still preferred for analyses run in parallel by chromosome. Step 3
ensures that each variant has a unique integer ID across the genome,
so the variant.id field in the per-chromosome files and the combined
file are consistent.


## Relatedness and Population structure

1. [KING-robust](http://www.ncbi.nlm.nih.gov/pubmed/20926424) to get
initial kinship estimates  
    `king.py`
    1. `ld_pruning.R`
    2. `combine_variants.R`
    3. `ibd_king.R`
2. [PC-AiR](http://www.ncbi.nlm.nih.gov/pubmed/25810074) to select an
informative set of unrelated samples, do PCA on unrelated, project
into relatives  
    `pcair.py`
    1. `find_unrelated.R`
    2. `ld_pruning.R`
    3. `combine_variants.R`
    4. `pca_byrel.R`
3. [PC-Relate](http://www.ncbi.nlm.nih.gov/pubmed/26748516) to
estimate kinship coefficients adjusted for population structure and
admixture using PCs  
    `pcrelate.py`
    1. `pcrelate.R`
4. Repeat steps 2-3, using new kinship estimates for PC-AiR


## Association testing

### Single-variant

1. `null_model.R`
2. `assoc_single.R`

### Aggregate

1. `null_model.R`
2. `assoc_window.R`
