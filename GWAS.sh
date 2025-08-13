#!/bin/bash

rm -f GWAS_outputs/* # quick clear of any previous runs

set -euo pipefail # prevents errors from being masked

cd GWAS_outputs

# ====== Defining ======
INFO="../data/20130606_g1k_3202_samples_ped_population.txt"   # 1000G metadata
CHRX="../data/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz" # ChrX
CHR="../data/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" # Chr22
SETS=("chrX" "chr22")  # Easy iteration of names (used in for loops)


# ====== FIXING Psam to contain correct info ======
#  ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
# | A text file which usually has at least one header line, where only the last header line starts with '#FID' or '#IID'. This final header line specifies the columns in the .psam file; the following intermediate column headers are recognized:
# |     IID (individual ID; required)
# |     SID (source ID, when there are multiple samples for the same individual)
# |     PAT (individual ID of father, '0' if unknown)
# |     MAT (individual ID of mother, '0' if unknown)
# |     SEX ('1' = male, '2' = female, 'NA'/'0' = unknown)
# |_______________________________________________________


# vv Setting up parents.txt and sex.txt to run --update for parents and sex. Best practice, safer than simply merging and dropping.
# AWK dataprocessing script to build parents.txt (FID IID PAT MAT) from INFO. Tab seperated
awk 'BEGIN{OFS="\t"}
  # skip header row
  NR==1 { next }
  {
    pat = ($3=="" || $3=="NA" || $3==".") ? 0 : $3
    mat = ($4=="" || $4=="NA" || $4==".") ? 0 : $4
    print $1, $2, pat, mat
  }' "${INFO}" > parents.txt  # checking INFO, set NA/./empty column 3 n 4 to 0. Must be in FID IID F M format.

# Same logic as above but to build sex.txt (FID IID SEX) from INFO
awk 'BEGIN{OFS="\t"} 
  NR==1 { next } 
  {
    sx = $5
    if (sx=="M" || sx=="m" || sx=="male" || sx=="MALE") sx = 1
    else if (sx=="F" || sx=="f" || sx=="female" || sx=="FEMALE") sx = 2
    else if (sx!="1" && sx!="2") sx = 0
    print $1, $2, sx
  }' "${INFO}" > sex.txt  # map sex to 1/2/0; in FID IID SEX format.

echo "sex head:\n"; head sex.txt
echo "parent head:\n"; head parent.txt #check column names and update IID and FID in CHRX

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

# checking: verify .psam exists for each set
for S in "${SETS[@]}"; do
  test -s "${S}.psam" || { echo "missing ${S}.psam"; exit 1; } # tests to see if chrX and chr22 created .psam correctly.
done  # guard before updates

echo "\nAll .psams generated!\n"

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

echo "\nparents update processed!\n"

# insert SID (Subject ID) column equal to IID (header: #FID IID SID PAT MAT SEX)
for S in "${SETS[@]}"; do
  awk 'BEGIN{OFS="\t"} 
    NR==1 {
      print "#FID", "IID", "SID", "PAT", "MAT", "SEX"
      next
    }
    {
      print $1, $2, $2, $3, $4, $5
    }' "${S}.psam" > "${S}.psam.tmp"  # add SID after IID
  
  mv "${S}.psam.tmp" "${S}.psam"  # finalize
done

echo "\nSID coloumn processed\n"

# Safety Checks.
wc -l chrX.psam chr22.psam  # row counts

diff \
  <(cut -f2 chrX.psam | tail -n +2) \ 
  <(cut -f2 chr22.psam | tail -n +2) \
  || echo "warning: sample order differs"  # tells if sample order netween ChrX and Chr22 are the same or not

for S in "${SETS[@]}"; do # for loop finding which files have which rows where IID != SID
  echo "Checking ${S}.psam:"
  awk 'NR>1 && $2!=$3 {print "  Mismatch row " NR ":", $2, "!=", $3}' "${S}.psam"
done