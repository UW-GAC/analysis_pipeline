# Array discordance pipeline for TOPMed

## Prepare array data

The following assumes that the array data is dbGaP fingerprint files with 10,000 SNPs.

1. Create directory for this study under `/projects/topmed/qc/compare_to_arrays/array_data/`. This is where we store results associated with the array data and not with any particular freeze. 
2. Convert tped to bed with `tped2bed.sh <file>`. bed file will be written in same directory as tped (under `/projects/topmed/downloaded_data/prior_array_data/`)
3. If necessary, liftover the .bim file to match the sequence build with `liftover_bim.R`.
4. Convert bed to gds with `qsub -N bed2gds runRscript.sh bed2gds.R <config>`. An AnnotatedDataFrame with sample annotation will be created also.

If we have an alternate array source, modify the above code so as to end up with a GDS file in the same build as the sequence data, and an AnnotatedDataFrame with columns "sample.id" and "subject.id". "sample.id" should be unique, and "subject.id" should match the "submitted_subject_id" in the TOPMed sample annotation.

## Compare array to sequence

5. Create directory for this study and the specific freeze to be checked under `/projects/topmed/qc/compare_to_arrays/<freeze>/`.
6. Subset sequencing GDS to overlapping variants and merge chromosomes: `array_subset.py <config>`

    config parameter | default value | description
    --- | --- | ---
    `out_prefix` | | Prefix for files created by this script.
    `array_gds_file` | | Path to array GDS file. 
    `seq_gds_file` | | Path to sequencing GDS file. 
	`seq_annot_file` | | Path to sequencing sample annotation file.
    `subset_gds_file` | | Path to output file.
    `study` | | Study name 
    `maf_threshold` | `0.05` | Minimum MAF for sequence variants to include
    `missing_threshold` | `0.01` | Maximum missing call rate for sequence variants to include

7. If we have an alternate array source, check the overlap between fingerprints and the subset file with `array_fingerprint_overlap.R`. If the number of overlapping variants is much less than 10,000, the code will supplement with a randomly selected set of variants that overlap between the sequence and array data, then combine these variants and the fingerprint variants in a GRanges object.
8. In the following steps, use a config file that sets `granges_include_file` to `dbgap_fp_ranges.hg38.RData`, or to the alternate include file defined in the previous step.
9. Run duplicate discordance with `array_disc.py -n N <config>` where N is the number of sample blocks to run in parallel.

    config parameter | default value | description
    --- | --- | ---
    `out_prefix` | | Prefix for files created by this script.
    `array_gds_file` | | Path to array GDS file.
	`array_annot_file` | | Path to array sample annotation file.
    `seq_gds_file` | | Path to sequencing GDS file (same as `subset_gds_file` above).
	`seq_annot_file` | | Path to sequencing sample annotation file.
    `study` | | Study name
	`granges_include_file` | `NA` | RData file with GRanges defining subset of variants.
	`sample_include_file` | `NA` | RData file with sample.id to include, if different from all samples in `study`.
	`variant_include_file` | `NA` | RData file with variants to include (result will be the intersection with `granges_include_file`)

10. If any samples were discordant, save an RData file with a vector of those sample ids. Check them against all other array samples with `qsub -t 1-N -N match_samples runRscript.sh -s array_match_sample.R <config>` where N is the number of samples to match.
