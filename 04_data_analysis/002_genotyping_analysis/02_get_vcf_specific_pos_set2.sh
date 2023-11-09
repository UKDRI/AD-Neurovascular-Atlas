#!/bin/bash

#SBATCH -p c_highmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=geneotyping_bbb_set1
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=20
#SBATCH --array=1-16%8
##### #SBATCH --mem-per-cpu=8000 # memory limit per core
#SBATCH --mem=100GB # memory limit per compute node for the job
#SBATCH --time=1-00:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/Hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/Hawk_output/%x_err_%A_%a_%J.txt

#### SBATCH --mem=740000 # memory limit per compute node for the job

echo "*****************************************************************"
echo "All jobs in this array have:"
echo "- SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID}"
echo "- SLURM_ARRAY_TASK_COUNT: ${SLURM_ARRAY_TASK_COUNT}"
echo "- SLURM_ARRAY_TASK_MIN: ${SLURM_ARRAY_TASK_MIN}"
echo "- SLURM_ARRAY_TASK_MAX: ${SLURM_ARRAY_TASK_MAX}"
echo "This job in the array has:"
echo "- SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "- SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Run on host: "`hostname`
echo "Number of threads (nproc): "`nproc`
echo "Total memory in GB: "`free -g | grep -oP '\d+' | sed -n 1p`
echo "Used memory in GB: "`free -g | grep -oP '\d+' | sed -n 2p`
echo "Free memory in GB: "`free -g | grep -oP '\d+' | sed -n 3p`
echo "Username: "`whoami`
echo "Started at: "`date`
echo -e "*****************************************************************\n"

module purge
module load samtools/1.10
module load vcftools/0.1.16
module load bcftools/1.10.2

# Pass the set number as an argument when running the script
# SET_NUM=$1

## BAM files from Cell Ranger count
BAM_DIR="03_data/990_processed_data/001_snrnaseq/04_cellranger_count/02_set2/"

## reference fasta file provided by Cell Ranger, part of reference download from website
CR_REF="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/03_cellranger_reference/GRCh38/fasta/genome.fa"

## sample IDs - correspond to names of Cell Ranger output folders
SAMPLE_ID_FILE="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/90_sample_info/samples_set2_rename.txt"

## output VCF files
OUTPUT_DIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/992_genotyping_data/01_vcf_specific_position/02_set2"

## 24 main chromosomes: as specified in Cell Ranger BAM files: chr1, chr2 ... chr22, chrX, chrY
MAIN_CHR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/992_genotyping_data/main_chromosomes_grch38_Cell_Ranger.txt"

## note, SNP positions on NeuroChip are annotated for GRCh37	
## used rtracklayer in R to liftover positions to hg38 reference
## from liftover_neurochip.Rmd
POS_FILE="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/992_genotyping_data/NeuroChip_v1_positions_liftover_GRCh38.txt"

#-----------------------------------------------------------------------

## samtools mpileup help

## Note that using "samtools mpileup" to generate BCF or VCF files is now
## deprecated.  To output these formats, please use "bcftools mpileup" instead.

# - --skip-indels still works but deprecated
# - option is still listed in bcftools mpileup

#-----------------------------------------------------------------------

mkdir -p $OUTPUT_DIR

N=${SLURM_ARRAY_TASK_ID}

SAMPLE_ID=$(cat $SAMPLE_ID_FILE | tail -n+${N} | head -1)

BAM_FILE=$BAM_DIR""$SAMPLE_ID"/outs/possorted_genome_bam.bam"
OUT_TMP=$OUTPUT_DIR""$SAMPLE_ID"_vcf_filenames.tmp"
OUT_IDX=$OUTPUT_DIR""$SAMPLE_ID
	
echo "Sample: "$SAMPLE_ID
echo "BAM file: "$BAM_FILE
echo "OUT idx: "$OUT_IDX

##  1. variant calling
# mpileup is rather slow, still not multi-threaded
# see request here: https://github.com/samtools/samtools/issues/480
# see also: https://github.com/brentp/bigly/issues/4
# but we are only interested in SNP array positions for genotype comparison
# hence perform mpileup only for those array positions
# faster with xargs parallel option, call mpileup for each main chromosome in parallel
cut -f1 $MAIN_CHR | uniq | xargs -I {} -n 1 -P 24 sh -c "samtools mpileup -uf $CR_REF --skip-indels -r {} -l $POS_FILE $BAM_FILE | bcftools call -mv > ${OUT_IDX}_var_raw_{}.vcf"

## 2. merge 24 individual VCF files into one VCF file per sample

# get header
grep "^#" ${OUT_IDX}"_var_raw_chr1.vcf" > ${OUT_IDX}"_var_raw.vcf"

cut -f1 $MAIN_CHR | uniq > ${OUT_TMP}

cat ${OUT_TMP} | while read TAG; do	
	FILE=${OUT_IDX}"_var_raw_"$TAG".vcf"	
	grep -v '^#' $FILE >> ${OUT_IDX}"_var_raw.vcf"
	rm $FILE
done

wait
rm $OUT_TMP

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
