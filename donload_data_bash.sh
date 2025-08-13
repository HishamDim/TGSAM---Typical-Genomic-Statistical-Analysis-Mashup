#!/bin/bash
cd data

dir="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/"
prefix="1kGP_high_coverage_Illumina.chr"
suffix=".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"


# === download for all chromosomes, (a lot of data so I'll just use chr22) ===
# for chr in {1..22} 
# do
#     wget "${dir}""${prefix}""${chr}""${suffix}".tbi "${dir}""${prefix}""${chr}""${suffix}"
# done
# # merge VCF files into a large file for all chormosomes
# bcftools concat chr1.vcf.gz chr2.vcf.gz ... chr22.vcf.gz -Oz -o merged.vcf.gz

# # index newly made file
# bcftools index merged.vcf.gz

# NOTE: THE CODE USES THE FOLLOWING ONLY, IF USING ALL CHR, ALL FILES WILL HAVE TO BE UPDATED ACCORDINGLY TO USE THE MERGED CHR MEGAFILE

# pop, family, and sex data (seperate in 1000genomes, this is atypical)
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

# this is X chromosome
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz

# just the 22nd chromosome for testing as its the smallest.
wget "${dir}""${prefix}""22""${suffix}".tbi "${dir}""${prefix}""22""${suffix}"