#!/bin/bash

# Following code uses plink2 to measure kinship of chr22 data

plink2 --vcf data/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --make-bed --out kin_temp_outputs/converted_output

plink2 --bfile kin_temp_outputs/converted_output --indep-pairwise 100 20 0.2 --out kin_temp_outputs/pruned_snps
plink2 --bfile kin_temp_outputs/converted_output --extract kin_temp_outputs/pruned_snps.prune.in --make-bed --out kin_temp_outputs/pruned_kin_data

plink2 --bfile kin_temp_outputs/pruned_kin_data --make-king-table --out kin_outputs/king_results