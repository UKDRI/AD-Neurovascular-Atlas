#!/bin/bash

#SBATCH -p c_compute_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=get_genotyping_comparison
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
##### #SBATCH --mem-per-cpu=8000 # memory limit per core
#SBATCH --mem=32GB # memory limit per compute node for the job
#SBATCH --time=0-10:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --mail-user Bernardo-HarringtonG@cardiff.ac.uk # email on fail
#SBATCH --mail-type END,FAIL

## 24 genotyped donors by NeuroChip
## HOW WERE THEY SPLIT? - looking at the VCFs, must be done via bcftools view
## 1 sample: 
##bcftools_viewCommand=view -s 1_NP002_2019 -o /Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/Genotyping/NeuroChip_24_results/Split_samples/1_NP002_2019.vcf /Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/Genotyping/NeuroChip_24_results/Webber_Neurochip_2021_chr.vcf; Date=Fri Jul 23 18:10:48 2021
INPUT_DIR_DNA="/Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/Genotyping/NeuroChip_24_results/Split_samples/"
INPUT_DIR_DNA_BATCH1="/gluster/dri02/rdscw/shared/webber/Endo_10X/Genotyping/NeuroChip_20_donors_set_1/Split_samples/"
INPUT_DIR_DNA_BATCH2="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/992_genotyping_data/PLINK_021123_0956/split_samples/"

## raw VCFs from Cell Ranger BAM files (lifted to hg19)
INPUT_DIR_RNA="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/992_genotyping_data/03_full_liftover_to_hg19/"

## store .gz VCF files and index
OUTPUT_DIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/992_genotyping_data/04_compare_samples_dna_rna_filtered/"

OUTPUT_DIR_DNA=$OUTPUT_DIR"VCF_DNA/"
OUTPUT_DIR_RNA=$OUTPUT_DIR"VCF_RNA/"

## results from pairwise comparison via gtcheck
## columns:
## sample1 sample2 n_sample1_passed n_sample2_passed n_discordance n_compared
OUT_FILE=$OUTPUT_DIR"gtcheck_40_samples_DP_20_800_qual_20.txt"

#BCFTOOLS="/Users/frankwessely/Desktop/Documents/Tools/bcftools-1.9/bcftools"
module purge
module load bcftools/1.10.2

#-----------------------------------------------------------------------

mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR_DNA
mkdir -p $OUTPUT_DIR_RNA

cd $OUTPUT_DIR

find $INPUT_DIR_RNA -maxdepth 2 -name '*_raw.vcf' > filenames_1.tmp

## 1. filter and index VCFs from RNA
cat filenames_1.tmp | while read VCF; do
	
	VCF_NAME=$(echo ${VCF##*/})
	VCF_GZ=$OUTPUT_DIR_RNA""$VCF_NAME".gz"
	
	echo $VCF
	echo $VCF_GZ
	echo " "
	
	## 1. filter and then compress VCF file
	## filter here somewhat arbitrary, could be optimised
	bcftools filter -s LowQual -e '%QUAL<20 || INFO/DP>800 || INFO/DP<20' $VCF | bcftools view -f PASS -O z -o $VCF_GZ
	
	## 2. index VCF file
	bcftools index $VCF_GZ
	
done

wait
rm filenames_1.tmp

#-----------------------------------------------------------------------

find $INPUT_DIR_DNA_BATCH1 -maxdepth 1 -name '*.vcf' > filenames_1.tmp
find $INPUT_DIR_DNA_BATCH2 -maxdepth 1 -name '*.vcf' >> filenames_1.tmp

## 2. index VCFs from DNA
cat filenames_1.tmp | while read VCF; do
	
	VCF_NAME=$(echo ${VCF##*/})
	VCF_GZ=$OUTPUT_DIR_DNA""$VCF_NAME".gz"
	
	echo $VCF
	echo $VCF_GZ
	echo " "
	
	## 1. compress VCF file
	bcftools view $VCF -O z -o $VCF_GZ

	## 2. index VCF file
	bcftools index $VCF_GZ
	
done

wait
rm filenames_1.tmp

#-----------------------------------------------------------------------

## 3. Compare VCFs via gtcheck

find $OUTPUT_DIR_DNA -maxdepth 1 -name '*.vcf.gz' > filenames_1.tmp
find $OUTPUT_DIR_RNA -maxdepth 1 -name '*.vcf.gz' > filenames_2.tmp

## new results file
echo -n > $OUT_FILE

cat filenames_1.tmp | while read VCF1; 
do
	cat filenames_2.tmp | while read VCF2; 
	do
	
		VCF_NAME_1=$(echo ${VCF1##*/} | sed 's/.vcf.gz//')
		VCF_NAME_2=$(echo ${VCF2##*/} | sed 's/.vcf.gz//')
		
		## pair of samples
		echo -e -n $VCF_NAME_1"\t"$VCF_NAME_2"\t" >> $OUT_FILE
		
		## number of records in each VCF (after above filtering)
		echo -e -n $(bcftools index --nrecords $VCF1)"\t"$(bcftools index --nrecords $VCF2)"\t" >> $OUT_FILE
		
		## gtcheck
		bcftools gtcheck -G 1 -g $VCF1 $VCF2 | grep '^CN' | cut -f2,4 >> $OUT_FILE
		
	done
done

wait
rm filenames_1.tmp
rm filenames_2.tmp
