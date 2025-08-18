#!/bin/bash

rm -f GWAS_outputs/* # quick clear of any previous runs

set -euo pipefail # prevents errors from being masked

cd GWAS_outputs

# TOC
# 1. Data preparing (Updating sex and putting in correct format)
# 2. Genotype Filtering
# 3. Sex Check
# 4. Remove Relatedness



# ====== Defining ======
INFO="../data/20130606_g1k_3202_samples_ped_population.txt"   # 1000G metadata
CHRX="../data/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz" # ChrX
CHR="../data/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" # Chr22
SETS=("chrX" "chr22")  # Easy iteration of names (used in for loops)


# ====== FIXING Psam to contain correct info ======
#  ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# | A text file which usually has at least one header line, where only the last header line starts with '#FID' or '#IID'. This final header line specifies the columns in the .psam file; the following intermediate column headers are recognized:
# |     #FID must either be the first column, or absent. If it's absent, all FID values are now assumed to be '0'.
# |     *#IID (individual ID; required)
# |     #SID (source ID, when there are multiple samples for the same individual)
# |     *PAT (individual ID of father, '0' if unknown)
# |     *MAT (individual ID of mother, '0' if unknown)
# |     *SEX ('1' = male, '2' = female, 'NA'/'0' = unknown)
# | * denotes required. Order is required.
# |_______________________________________________________


# vv Setting up parents.txt and sex.txt to run --update for parents and sex. Best practice, note that the formats are REQUIRED by PLINK2
# AWK dataprocessing script to build parents.txt (IID FID PAT MAT) from INFO. Tab seperated.
awk 'BEGIN{ OFS="\t" }
# skip header row
NR==1 { print "#IID", "#FID", "PAT", "MAT"; next}
{
  pat = ($3=="" || $3=="NA" || $3==".") ? 0 : $3
  mat = ($4=="" || $4=="NA" || $4==".") ? 0 : $4
  print $2, $1, pat, mat
}' "${INFO}" > parents.txt  # checking INFO, set NA/./empty column 3 n 4 to 0. Must be in IID FID M F format according to original info data.

# Same logic as above but to build sex.txt (IID FID SEX) from INFO
awk 'BEGIN{OFS="\t"}
NR==1 { print "#IID","#FID","SEX"; next }
{
  print $2, $1, substr($5,1,1)
}' "${INFO}" > sex.txt
# map sex to 1/2/0; in IID FID SEX format.

echo "sex head:"; head sex.txt
echo "parent head:"; head parents.txt

# (====== MAKING Pgen Pvar and Psam for X and 22 Chr. ======)
plink2 --vcf "${CHRX}" \
  --update-sex sex.txt \
  --split-par b37 \
  --make-pgen \
  --out chrX
# Note update sex was run in file formation for ChrX, not like Chr22 below as plink2 requires sex information.
# split-par b37 denotes the build of the dataset. Chromosome X has two pseudoautosomal regions (PAR1 and PAR2) that are homologous to regions on chromosome Y and behave like autosomal regions during recombination. PLINK2 requires the --split-par option to properly handle these regions when they are present in the VCF file, as it splits PAR and non-PAR variants into separate categories (e.g., assigning PAR variants to a special chromosome code like PAR1 or PAR2).

plink2 --vcf "${CHR}" \
  --update-sex sex.txt \
  --make-pgen \
  --out chr22
# Same logic as before, sex update not required here, but added for compactness isntead of being done seperatly later.

head chrX.psam
head chr22.psam

# checking: verify .psam exists for each set
for S in "${SETS[@]}"; do
  test -s "${S}.psam" || { echo "missing ${S}.psam"; exit 1; } # tests to see if chrX and chr22 created .psam correctly.
done  # guard before updates

echo "All .psams generated!"

# apply PLINK2s default updates and rewrite pgen/pvar/psam for PARENTS only (sex already done)
for S in "${SETS[@]}"; do
  plink2 -pfile "${S}" \
    --update-parents parents.txt \
    --make-pgen \
    --out "${S}.tmp"  # fill PAT/MAT and persist, .tmp to prevent corruption of og data during update.
  
  mv "${S}.tmp.pgen" "${S}.pgen"
  mv "${S}.tmp.pvar" "${S}.pvar"
  mv "${S}.tmp.psam" "${S}.psam"  # replace in place
done

echo "parents update processed!"
echo "chr22:"; head chr22.psam
echo "chrX:"; head chrX.psam

# NOTE: #FID and #SID are not used in this analysis as they are optional, but if used they must follow the order above.

# Safety Checks.
wc -l chrX.psam chr22.psam  # row counts



# ====== 2) Basic genotype filtering ======

QC_MAF=0.01 # variant-level: keep MAF >= 1% (minor allele frequency; too few carriers in the dataset to reliably detect associations)
QC_GENO=0.02 # variant-level: keep variants with missingness <= 2%. | FOR REMOVING SNPs & VARIANTS MISSING TOO MANY PEOPLE
QC_MIND=0.02 # sample-level: remove samples with missingness > 2% | FOR REMOVING PEOPLE W/ TOO MANY MISSING GENOTYPES

# calculates missingness statistics
for S in "${SETS[@]}"; do
  plink2 \
    --pfile "${S}" \
    --missing sample-only \
    --out "${S}.miss"  # writes ${S}.miss.smiss with F_MISS per sample
done

# Assigns QC_MIND to var t and uses it to compare values of missingness and only transcribe those that exceed t into mind_fail.txt
# for loop identifies columns (dynamic and adaptive, does not require changing in different tables, unlike prior)
awk -v t="$QC_MIND" '
FNR==1 {
  if (NR==1) {
    for (i=1; i<=NF; i++){
      if ($i == "F_MISS") c=i
      if ($i == "#IID") ci=i
    }
    print "#IID"
  }
  next
}
($c+0) > t { print $ci }
' "${S}".miss.smiss \
  | sort -u >> mind_fail.txt # run awk on both sample missingness files n sort IIDs and remove duplicates, append to mind_fail.txt
# didn't use for loop because appending not making two seperate files. (Even though removing dupes, for best practice)

echo "made mind_fail.txt"

# apply filters per dataset
for S in "${SETS[@]}"; do
  plink2 \
    --pfile "${S}" \
    --remove mind_fail.txt \
    --maf "${QC_MAF}" \
    --geno "${QC_GENO}" \
    --make-pgen \
    --out "${S}.qc"  # outputs: ${S}.qc.{pgen,pvar,psam}
done

echo -e "\nDONE ALL BASIC GENOTYPE FILTERING"

# quick counts
wc -l mind_fail.txt  # number of removed samples
for S in "${SETS[@]}"; do
  echo "${S}.qc variants:"; awk 'END{print NR-1}' "${S}.qc.pvar"  # pvar lines minus header
done





# ====== 3) Sex check (X-chromosome F-statistic against reported SEX)======
echo -e "Beginning Sex Check"
# cd GWAS_outputs

plink2 \
  --pfile chrX.qc \
  --check-sex \
  --out sexchk  # writes sexchk.sexcheck

awk 'NR>1 {print $5}' sexchk.sexcheck | sort -n > F_values.txt  # look at F values
# Does a basic sex check and writes the F_values to a text file to be examined in "Visualizations". This allows for customization of F_value thresholds.
# exit

rm -f sexchk.sexcheck # removes old sexcheck
rm -f sexchk.log # removes old sexcheck

SEX_LO=0.5 # lower F cutoff for females
SEX_HI=0.8 # upper F cutoff for males

# Lets re-run it with proper threshold values
plink2 \
  --pfile chrX.qc \
  --check-sex max-female-xf=${SEX_LO} min-male-xf=${SEX_HI} \
  --out sexchk_lo${SEX_LO}_hi${SEX_HI} # writes .sexcheck

awk 'NR==1 {
  for(i=1; i<=NF; i++) {h[$i]=i}
  print $h["#IID"],$h["PEDSEX"],$h["SNPSEX"],$h["F"],$h["STATUS"]
  next
}
$h["STATUS"]!="OK" {
  print $h["#IID"],$h["PEDSEX"],$h["SNPSEX"],$h["F"],$h["STATUS"]
}' OFS="\t" sexchk_lo${SEX_LO}_hi${SEX_HI}.sexcheck > sex_mismatch.tsv  # audit table: IID, reported SEX, SNP-inferred sex, F, status

awk 'NR==1{print "#IID"; next} {print $1}' sex_mismatch.tsv > sex_fail.remove  # IID-only remove list for downstream filters

wc -l sex_mismatch.tsv                          # count flagged samples (includes header)
head -n 10 sex_mismatch.tsv                     # quick peek at top rows

for S in "${SETS[@]}"; do
  plink2 \
    --pfile "${S}.qc" \
    --remove sex_fail.remove \
    --make-pgen \
    --out "${S}.nosexmismatch"
done

# ====== 4) Removing Kin ======