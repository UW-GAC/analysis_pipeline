# TOPMed analysis pipeline

## Setup

1. Install R packages and dependencies from Bioconductor
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("SeqVarTools", "SNPRelate", "GENESIS"))
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

## Conversion to GDS

1. vcf2gds.R
2. merge_gds.R
3. unique_variant_ids.R

Step 1 converts VCF files (one per chromosome) into GDS files,
discarding non-genotype FORMAT fields. Step 2 combines these files
into a single GDS file, which is needed for whole-genome analyses such
as relatedness and population structure. The single-chromosome files
are still preferred for analyses run in parallel by chromosome. Step 3
ensures that each variant has a unique integer ID across the genome,
so the variant.id field in the per-chromosome files and the combined
file are consistent.

## Relatedness and Population structure

1. ld_pruning.R
2. ibd_king.R
3. pcair.R (or pca_snprelate.R)
4. pcrelate.R

## Association testing

### Single-variant

1. null_model.R
2. assoc_single.R

### Aggregate

1. null_model.R
2. assoc_window.R
