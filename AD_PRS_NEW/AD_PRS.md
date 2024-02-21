

## Ref

[PRS Pipeline by Jennifer Collister, Xiaonan Liu](https://2cjenn.github.io/PRS_Pipeline/)


## Load modules

```bash
# bgen module 
module load bgen
bgenix -help
cat-bgen -help

# plink module
module load plink
plink2 --version

# SQLite module
module load SQLite
sqlite3 --help
```

## Extract required SNPs

```BASH
# SNP list
cat Betas_ADPGS_final.csv | awk -F, 'NR>1{ print $1 }' > rsidlist.txt
cat Betas_ADPGS_final.csv | awk -F, 'NR>1{ print sprintf("%02d", $2)":"$3"-"$3 }' > chrposlist.txt
# Process bgen file for one chrom
i=17
bgenix -g /data/LEPSUKB/UKB_imputation_from_genotype/ukb22828_c${i}_b0_v3.bgen -incl-rsids rsidlist.txt -incl-range chrposlist.txt > bgen/chr_${i}.bgen
# Process bgen file for multiple chrom
#for i in {19..22}; do
for i in {1..22}; do
    echo "############################################"
    echo "############################################"
    echo "Processing chr${i}: $(date)"
    bgenix \
        -g /data/LEPSUKB/UKB_imputation_from_genotype/ukb22828_c${i}_b0_v3.bgen \
        -incl-rsids rsidlist.txt \
        -incl-range chrposlist.txt \
        > bgen/chr_${i}.bgen
done

# Combine the .bgen files for each chromosome into one
## For unknow reasons, cat-bgen cannot expand bash variables, $bg_file is interpreted
## literally as $bg_file string. So the following strategy uses echo to 
## expand the variable first. 
## Generate bgen file list, replace newline with space
bg_file=$(ls -v bgen/*.bgen | paste -s -d ' ')
echo $bg_file
## Visually check the command
echo "cat-bgen -g  $bg_file -og multiple_chr.bgen -clobber"  
## Run it
echo "cat-bgen -g  $bg_file -og bgen/multiple_chr.bgen -clobber" | bash 

# Generate index
bgenix -g bgen/multiple_chr.bgen -index -clobber
# Check bgen file 
bgenix -g bgen/multiple_chr.bgen -list | less -SN
bgenix -g bgen/multiple_chr.bgen -vcf | less -SN
```

## SNP alignment and data conversion

```bash
cd bgen
# Import the betas into the sqlite database as a table called Betas
## Drop Betas table from multiple_chr.bgen.bgi database
sqlite3 multiple_chr.bgen.bgi "DROP TABLE IF EXISTS Betas;"
## Import csv file into database as Betas table 
sqlite3 -separator "," multiple_chr.bgen.bgi ".import ../Betas_ADPGS_final.csv Betas"
## Drop Joined table from multiple_chr.bgen.bgi database
sqlite3 multiple_chr.bgen.bgi "DROP TABLE IF EXISTS Joined;"
## Some basic commands
sqlite3 multiple_chr.bgen.bgi ".help"
sqlite3 multiple_chr.bgen.bgi ".show"
sqlite3 multiple_chr.bgen.bgi ".tables"
sqlite3 multiple_chr.bgen.bgi ".schema Betas"
sqlite3 multiple_chr.bgen.bgi ".schema Variant"

# And inner join it to the index table (Variants), making a new table (Joined)
# By joining on alleles as well as chromosome and position 
# we can ensure only the relevant alleles from any multi-allelic SNPs are retained
sqlite3 -header -csv multiple_chr.bgen.bgi \
"CREATE TABLE Joined AS 
  SELECT Variant.*, Betas.chr_name, Betas.Beta FROM Variant INNER JOIN Betas 
    ON Variant.chromosome = printf('%02d', Betas.chr_name) 
    AND Variant.position = Betas.chr_position 
    AND Variant.allele1 = Betas.noneffect_allele 
    AND Variant.allele2 = Betas.effect_allele 
  UNION 
  SELECT Variant.*, Betas.chr_name, -Betas.Beta FROM Variant INNER JOIN Betas 
    ON Variant.chromosome = printf('%02d', Betas.chr_name) 
    AND Variant.position = Betas.chr_position 
    AND Variant.allele1 = Betas.effect_allele 
    AND  Variant.allele2 = Betas.noneffect_allele;"

# Filter the .bgen file to include only the alleles specified in the Betas for each SNP 
bgenix -g multiple_chr.bgen -table Joined  > single_allelic.bgen
# And produce an index file for the new .bgen
bgenix -g single_allelic.bgen -index



# File conversion
# Load the data from .bgen and .sample files, and convert it to PLINK 2 format.
#--sample ukbA_imp_chrN_v3_sP.sample \
plink2 --bgen single_allelic.bgen ref-first \
--hard-call-threshold 0.1 \
--sample /data/LEPSUKB/UKB_imputation_from_genotype/ukb22828_c1_b0_v3.sample \
--memory 15000 \
--set-all-var-ids @:#_\$r_\$a \
--freq \
--make-pgen \
--out raw 

```

## QC 

```bash
# Ambiguous SNPs
## Talindromic SNPs (ie SNPs where both alleles are (A, T) or (C, G)) which have effect allele frequencies close enough to 0.5 that we cannot unambiguously distinguish between the effect and non-effect alleles.
awk '/^[^#]/ { if( $5>0.4 && $5<0.6 && ( ($3=="A" && $4=="T") || ($4=="T" && $3=="A") || ($3=="C" && $4=="G") || ($4=="G" && $3=="C") ) ) { print $0 }}' \
raw.afreq > exclrsIDs_ambiguous.txt

# Imputation information
for i in {1..22}; do
#for i in {20..22}; do
  awk -v chr=$i 'BEGIN {FS="\t"; OFS="\t"} { print chr,$0,chr":"$3"_"$4"_"$5 }' /data/LEPSUKB/UKB_imputation_from_genotype/ukb22828_c${i}_b0_v3.mfi.txt
done > ukb_mfi_all_v3.tsv

# Filter SNPs
plink2 --pfile raw \
--memory 15000 \
--exclude exclrsIDs_ambiguous.txt \
--extract-col-cond ukb_mfi_all_v3.tsv 9 10 --extract-col-cond-min 0.4 \
--maf 0.01 \
--write-snplist \
--make-pgen \
--out snpQC 

# Sample QC
plink2 --pfile raw \
--memory 15000 \
--extract snpQC.snplist \
--keep-fam ../usedinpca.txt \
--write-samples \
--out sampleQC 

```

## Calculate the PRS

```bash
# Prepare the score.txt file. 
# Rearrange the SNPs and betas from the Betas.csv file into the format the PLINK 2 --score command expects.
sqlite3 -separator " " -list multiple_chr.bgen.bgi \
"SELECT chr_name || ':' || position || '_' ||allele1 || '_' || allele2, allele2, Beta FROM Joined;" > score.txt

# Run plink2 --score command
plink2 --pfile raw \
--memory 15000 \
--extract snpQC.snplist \
--keep sampleQC.id \
--score score.txt no-mean-imputation \
--out PRS

```




```bash


```




