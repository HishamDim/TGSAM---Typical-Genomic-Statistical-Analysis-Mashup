#!/bin/bash

# Following code performs a PCA on Chr22

# Clears any previous runs
rm -f PCAoutputs/*
rm -f PCAtemp_outputs/*

# convert VCF into PLINK format
plink2 --vcf data/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz --make-pgen \
    --out PCAtemp_outputs/converted_output

# === Pruning to remove high correlated and prevent clustering ===
plink2 --pfile PCAtemp_outputs/converted_output --rm-dup --snps-only --maf 0.05 --indep-pairwise 100 20 0.2 \
    --out PCAtemp_outputs/pruned_SNPs_tokeep
plink2 --pfile PCAtemp_outputs/converted_output --extract PCAtemp_outputs/pruned_SNPs_tokeep.prune.in --king-cutoff 0.0884 \
    --make-pgen --out PCAtemp_outputs/output_pruned
# 0.0884 is first degree cousin, meaning it doesnt include anything MORE related than first cousin

# Running PCA
plink2 --pfile PCAtemp_outputs/output_pruned --pca 20 --out PCAoutputs/pca_results