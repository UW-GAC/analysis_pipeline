# TOPMed analysis pipeline

## Setup

1. Install R packages and dependencies from Bioconductor  
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("SeqVarTools", "SNPRelate", "GENESIS", "argparser", "dplyr", "tidyr", "ggplot2", "GGally"))
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

Some scripts can be run in parallel by chromosome. For these scripts, the chromosome number is the second argument. If running in parallel, include a space in file names in the config file where chromosome should be inserted, e.g.,
```
gds_file "1KG_phase3_subset_chr .gds"
```

Nearly all scripts require a GDS file in SeqArray format. Phenotype files should be an AnnotatedDataFrame saved in an RData file. See `?AnnotatedDataFrame` or the SeqVarTools documentation for details. Example files are provided in `testdata`.

Python scripts are provided to run multi-step analyses on a compute cluster or cloud environment. The `submitJob` function in `TopmedPipeline.py` is written for a Sun Grid Engine (SGE) cluster, but may be modified for other environments. These scripts require a config argument `out_prefix` in addition to the arguments for each R script called. Some input and output file name parameters are overridden by the scripts in order to link jobs together. Example config files are in `testdata`.


## Conversion to GDS

`vcf2gds.py`

1. `vcf2gds.R`
2. `merge_gds.R`
3. `unique_variant_ids.R`

config parameter | default value | description
--- | --- | ---
`out_prefix` | | Prefix for files created by this script.
`vcf_file` | | Input VCF file. Include a space to insert chromosome.
`gds_file` | | Output GDS file. Include a space to insert chromosome.
`merged_gds_file` | | Merged genotype-only GDS file containing all chromosomes.
`format` | `GT` | FORMAT fields from the VCF to convert to GDS. Default is genotypes only.

Step 1 converts VCF files (one per chromosome) into GDS files, discarding non-genotype FORMAT fields. Step 2 combines these files into a single GDS file, which is needed for whole-genome analyses such as relatedness and population structure. The single-chromosome files are still preferred for analyses run in parallel by chromosome. Step 3 ensures that each variant has a unique integer ID across the genome, so the variant.id field in the per-chromosome files and the combined file are consistent.


## Relatedness and Population structure

1. [KING-robust](http://www.ncbi.nlm.nih.gov/pubmed/20926424) to get initial kinship estimates

    `king.py`
    1. `ld_pruning.R`
    2. `combine_variants.R`
    3. `ibd_king.R`
	4. `kinship_plots.R`

    config parameter | default value | description
    --- | --- | ---
    `out_prefix` | | Prefix for files created by this script.
    `gds_file` | | GDS file with all chromosomes.
	`ld_r_threshold` | `0.32` | `r` threshold for LD pruning. Default is `r^2 = 0.1`.
	`ld_win_size` | `10` | Sliding window size in Mb for LD pruning.
	`maf_threshold` | `0.01` | Minimum MAF for variants used in LD pruning.
	`sample_include_file` | `NA` | RData file with vector of sample.id to include.
	`variant_include_file` | `NA` | RData file with vector of variant.id to consider for LD pruning. 
	`phenotype_file` | `NA` | RData file with AnnotatedDataFrame of phenotypes. Used for plotting kinship estimates separately by study.
	`study` | `NA` | Name of column in `phenotype_file` containing study variable.

2. [PC-AiR](http://www.ncbi.nlm.nih.gov/pubmed/25810074) to select an informative set of unrelated samples, do PCA on unrelated, project into relatives

    `pcair.py`
    1. `find_unrelated.R`
    2. `ld_pruning.R`
    3. `combine_variants.R`
    4. `pca_byrel.R`
	5. `pca_plots.R`
	6. `pca_corr.R`
	7. `pca_corr_plots.R`
	
    config parameter | default value | description
    --- | --- | ---
    `out_prefix` | | Prefix for files created by this script.
    `gds_file` | | GDS file with all chromosomes.
	`king_file` | | RData file with kinship coefficients created by `king.py`.
	`kinship_method` | `king` | Type of kinship estimates to use for finding unrelated set. Options are `king` or `pcrelate`.
	`kinship_threshold` | `0.04419417` | Minimum kinship estimate to use for assigning relatives (default is `2^(-9/2)` or 3rd degree relatives).
	`pcrelate_file` | `NA` | GDS file created by `pcrelate.py`. Only used if `kinship_method` is `pcrelate`.
	`sample_include_file` | `NA` | RData file with vector of sample.id to include.
	`ld_r_threshold` | `0.32` | `r` threshold for LD pruning. Default is `r^2 = 0.1`.
	`ld_win_size` | `10` | Sliding window size in Mb for LD pruning.
	`maf_threshold` | `0.01` | Minimum MAF for variants used in LD pruning.
	`variant_include_file` | `NA` | RData file with vector of variant.id to consider for LD pruning. 
	`n_pcs` | `20` | Number of PCs to return. 
	`n_pair` | `6` | Number of PCs in include in the pairs plot.
	`n_perpage` | `4` | Number of PC-variant correlation plots to stack in a single page. The number of png files generated will be `ceiling(n_pcs/n_perpage)`. 
	`thin` | `TRUE` | Logical for whether to thin points in the PC-variant correlation plots.
	`phenotype_file` | `NA` | RData file with AnnotatedDataFrame of phenotypes. Used for color-coding PCA plots by group.
	`group` | `NA` | Name of column in `phenotype_file` containing group variable.

3. [PC-Relate](http://www.ncbi.nlm.nih.gov/pubmed/26748516) to estimate kinship coefficients adjusted for population structure and admixture using PCs

    `pcrelate.py`
    1. `pcrelate.R`
	2. `kinship_plots.R`
	
    config parameter | default value | description
    --- | --- | ---
    `out_prefix` | | Prefix for files created by this script. 
    `gds_file` | | GDS file with all chromosomes.
	`pca_file` | | RData file with PCA results created by `pcair.py`. 
	`variant_include_file` | | RData file with LD pruned variant.id created by `pcair.py`.
	`n_pcs` | `3` | Number of PCs to use in adjusting for ancestry.
	`sample_block_size` | `10000` | Maximum number of samples to read in a single block. Adjust depending on computer memory and number of samples in the analysis. 
	`sample_include_file` | `NA` | RData file with vector of sample.id to include. 
	`phenotype_file` | `NA` | RData file with AnnotatedDataFrame of phenotypes. Used for plotting kinship estimates separately by study.
	`study` | `NA` | Name of column in `phenotype_file` containing study variable.
	
4. Repeat steps 2-3, using new kinship estimates for PC-AiR


## Association testing

An inverse-normal transform may be requested with `inverse_normal TRUE` in the config file. This is done by fitting the null model and rank-normalizing the marginal residuals. The normalized residuals are then used as the outcome in a null model with only PCs and kinship as covariates.

### Single-variant

`assoc.py single` 

1. `null_model.R`
2. `assoc_single.R`
3. `assoc_plots.R`

config parameter | default value | description
--- | --- | ---
`out_prefix` | | Prefix for files created by this script. 
`gds_file` | | GDS file. Include a space to insert chromosome.
`pca_file` | | RData file with PCA results created by `pcair.py`. 
`pcrelate_file` | | GDS file created by `pcrelate.py`.
`phenotype_file` | | RData file with AnnotatedDataFrame of phenotypes. 
`outcome` | | Name of column in `phenotype_file` containing outcome variable.
`binary` | `FALSE` | `TRUE` if `outcome` is a binary (case/control) variable; `FALSE` if `outcome` is a continuous variable.
`covars` | `NA` | Names of columns `phenotype_file` containing covariates, quoted and separated by spaces.
`inverse_normal` | `FALSE` | `TRUE` if an inverse-normal transform should be applied to the outcome variable.
`n_pcs` | `3` | Number of PCs to include as covariates.
`sample_include_file` | `NA` | RData file with vector of sample.id to include.
`mac_threshold` | `5` | Minimum minor allele count for variants to include in test. Use a higher threshold when outcome is binary.
`maf_threshold` | `0.001` | Minimum minor allele frequency for variants to include in test. Only used if `mac_threshold` is `NA`.
`pass_only` | `TRUE` | `TRUE` to select only variants with FILTER=PASS.
`test_type` | `score` | Type of test to perform. Options are `score` and `wald` if `binary` is `FALSE`, `score` only if `binary` is `TRUE`.
`variant_include_file` | `NA` | RData file with vector of variant.id to include.

### Aggregate

`assoc.py aggregate` 

1. `null_model.R`
2. `aggregate_list.R`
3. `assoc_aggregate.R`
4. `assoc_plots.R`

config parameter | default value | description
--- | --- | ---
`out_prefix` | | Prefix for files created by this script. 
`gds_file` | | GDS file. Include a space to insert chromosome.
`aggregate_type` | `allele` | Type of aggregate grouping. Options are to select variants by `allele` (unique variants) or `position` (regions of interest).
`variant_group_file` | | RData file with data frame defining aggregate groups. If `aggregate_type` is `allele`, columns should be `group_id`, `chromosome`, `position`, `ref`, `alt`. If `aggregate_type` is `position`, columns should be `group_id`, `chromosome`, `start`, `end`.
`pca_file` | | RData file with PCA results created by `pcair.py`. 
`pcrelate_file` | | GDS file created by `pcrelate.py`.
`phenotype_file` | | RData file with AnnotatedDataFrame of phenotypes. 
`outcome` | | Name of column in `phenotype_file` containing outcome variable.
`binary` | `FALSE` | `TRUE` if `outcome` is a binary (case/control) variable; `FALSE` if `outcome` is a continuous variable.
`covars` | `NA` | Names of columns `phenotype_file` containing covariates, quoted and separated by spaces. 
`inverse_normal` | `FALSE` | `TRUE` if an inverse-normal transform should be applied to the outcome variable.
`n_pcs` | `3` | Number of PCs to include as covariates.
`sample_include_file` | `NA` | RData file with vector of sample.id to include. 
`variant_include_file` | `NA` | RData file with vector of variant.id to include. Variants used will be the intersection of this set and variants defined by `variant_group_file`.
`alt_freq_range` | `"0 1"` | Range of alternate allele frequencies to consider, quoted and separated by spaces.
`test` | `burden` | Test to perform. Options are `burden` or `skat`.
`test_type` | `score` | Type of test to perform if `test` is `burden`. Options are `score` and `wald` if `binary` is `FALSE`, `score` only if `binary` is `TRUE`. 
`pval_skat` | `kuonen` | Method used to calculate p-values if `test` is `skat`. Options are `kuonen` (uses saddlepoint method), `davies` (uses numerical integration), and `liu` (uses a moment matching approximation). 
`rho` | `0` | A numeric value (or quoted, space-delimited list of numeric values) in [0,1] specifying the rho parameter when `test` is `skat`. `0` is a standard SKAT test, `1` is a score burden test, and multiple values is a SKAT-O test.
`weight_beta` | `"0.5 0.5"` | Parameters of the Beta distribution used to determine variant weights, quoted and space-delimited. `"0.5 0.5"` is proportional to the Madsen-Browning weights and `"1 25"` gives the Wu weights.

### Sliding window

`assoc.py window` 

1. `null_model.R`
2. `assoc_window.R`
3. `assoc_plots.R`

config parameter | default value | description
--- | --- | ---
`out_prefix` | | Prefix for files created by this script. 
`gds_file` | | GDS file. Include a space to insert chromosome.
`pca_file` | | RData file with PCA results created by `pcair.py`. 
`pcrelate_file` | | GDS file created by `pcrelate.py`.
`phenotype_file` | | RData file with AnnotatedDataFrame of phenotypes. 
`outcome` | | Name of column in `phenotype_file` containing outcome variable.
`binary` | `FALSE` | `TRUE` if `outcome` is a binary (case/control) variable; `FALSE` if `outcome` is a continuous variable.
`covars` | `NA` | Names of columns `phenotype_file` containing covariates, quoted and separated by spaces. 
`inverse_normal` | `FALSE` | `TRUE` if an inverse-normal transform should be applied to the outcome variable.
`n_pcs` | `3` | Number of PCs to include as covariates.
`sample_include_file` | `NA` | RData file with vector of sample.id to include. 
`variant_include_file` | `NA` | RData file with vector of variant.id to include. 
`alt_freq_range` | `"0 1"` | Range of alternate allele frequencies to consider, quoted and separated by spaces.
`test` | `burden` | Test to perform. Options are `burden` or `skat`. 
`test_type` | `score` | Type of test to perform if `test` is `burden`. Options are `score` and `wald` if `binary` is `FALSE`, `score` only if `binary` is `TRUE`. 
`pval_skat` | `kuonen` | Method used to calculate p-values if `test` is `skat`. Options are `kuonen` (uses saddlepoint method), `davies` (uses numerical integration), and `liu` (uses a moment matching approximation). 
`rho` | `0` | A numeric value (or quoted, space-delimited list of numeric values) in [0,1] specifying the rho parameter when `test` is `skat`. `0` is a standard SKAT test, `1` is a score burden test, and multiple values is a SKAT-O test.
`weight_beta` | `"0.5 0.5"` | Parameters of the Beta distribution used to determine variant weights, quoted and space-delimited. `"0.5 0.5"` is proportional to the Madsen-Browning weights and `"1 25"` gives the Wu weights.
`window_size` | `50` | Size of sliding window in kb. 
`window_step` | `20` | Step size of sliding window in kb.
