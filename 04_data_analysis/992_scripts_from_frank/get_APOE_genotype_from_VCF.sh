#!/bin/bash

## get APOE genotypes
## 
## prepare VCF to be used with apoe-genotyper.py
## available from: https://github.com/jjfarrell/apoe-genotyper
##
## NOTE: python script will use the last bead probe (row) reported in VCF file at the 2 APOE SNP positions
## script could be updated to check all or the subset of probes that should perform well
## apparently, some probes might not perform well:
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5534378/
## Given the importance of APOE, NeuroChip was designed so that rs7412 is genotyped by four separate probes 
## (three of which performed well: rs7412, seq-rs7412-B1, seq-rs7412-B3). 
## Similarly, rs429358 was genotyped by five separate bead probes (two of which performed well: seq-rs429358-T2, seq-rs429358-T3).
##
## the last probes in the VCF file: 
## seq-rs7412-T2 (actually NOT mentioned to perform well, but had the same calls as the other 3, so better check probes)
## seq-rs429358-T3 (should perform well)
##
## dependencies:
## pip install PyVCF
## pip install pysam

#-----------------------------------------------------------------------

## 24 genotyped donors by NeuroChip (4 are replicates to fill the chip, suffix '_R')
## based on GRCh37
INPUT_VCF="/Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/Genotyping/NeuroChip_24_results/Webber_Neurochip_2021.vcf"

## store compressed VCF, index and genotyping results
OUT_DIR="/Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/Genotyping/APOE_genotypes/"

## .gz out file
OUT_VCF="Webber_Neurochip_2021.vcf.gz"

## APOE genotyper python 3 script
APOE_GENOTYPER="/Users/frankwessely/Desktop/Documents/Tools/APOE_genotyper/apoe-genotyper.py"

## bcftools to generate compressed VCF with tabix index
BCFTOOLS="/Users/frankwessely/Desktop/Documents/Tools/bcftools-1.9/bcftools"

#-----------------------------------------------------------------------

mkdir -p $OUT_DIR

## 1. compress and index VCF
$BCFTOOLS view $INPUT_VCF -O z -o $OUT_DIR""$OUT_VCF
$BCFTOOLS index --tbi $OUT_DIR""$OUT_VCF

## 2. run APOE genotyper
python $APOE_GENOTYPER --vcf $OUT_DIR""$OUT_VCF --genome GRCh37 --project "endo_10X_24_samples" --out $OUT_DIR"apoe_24_neurochip"


#-----------------------------------------------------------------------
## check whether 2 or 4 probes give identical genotype calls

## grep -E '45411941|45412079' $INPUT_VCF

echo -e "\nNeuroChip probes check:"

echo -e 'seq-rs429358-T2\nseq-rs429358-T3\nrs7412\nseq-rs7412-B1\nseq-rs7412-B3\nseq-rs7412-T2' > dummy.txt
$BCFTOOLS filter -i 'ID=@dummy.txt' $INPUT_VCF > dummy.vcf

## check only 2 probes (out of 5)
P1=$(grep 'seq-rs429358-T2' dummy.vcf | cut -f 10-)
P2=$(grep 'seq-rs429358-T3' dummy.vcf | cut -f 10-)

if [[ "$P1" == "$P2" ]]
then
	echo "seq-rs429358-T2 and seq-rs429358-T3: identical"
else
	echo "seq-rs429358-T2 and seq-rs429358-T3: not identical"
fi

## check all 4 probes
P1=$(grep $'\trs7412\t' dummy.vcf | cut -f 10-)
P2=$(grep 'seq-rs7412-B1' dummy.vcf | cut -f 10-)
P3=$(grep 'seq-rs7412-B3' dummy.vcf | cut -f 10-)
P4=$(grep 'seq-rs7412-T2' dummy.vcf | cut -f 10-)

if [[ "$P1" == "$P2" && "$P1" == "$P3" && "$P1" == "$P4" ]]
then
	echo "4 rs7412 probes: identical"
else
	echo "4 rs7412 probes: not identical"
fi

rm dummy.txt
rm dummy.vcf

