#!/bin/bash

#SBATCH -p c_highmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=cellranger_count_set2
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-16%8
##### #SBATCH --mem-per-cpu=8000 # memory limit per core
#SBATCH --mem=340G # memory limit per compute node for the job
#SBATCH --time=3-00:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --mail-user Bernardo-HarringtonG@cardiff.ac.uk # email on fail
#SBATCH --mail-type END,FAIL

## Project root variables defined in the following
source 00_hpc_variables.sh

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


## FASTQ files, set 2 with 16 samples, data from 4 runs, sequenced in Oxford
## NOTE, V_19, only sequenced in 1 run (INPUT_DIR_1)
## NOTE, V_20, 1 run (INPUT_DIR_2) only with a single read, skip that run
## handled via if statement below
INPUT_DIR_1=$INPUT_DIR_b1s2"210707_A00711_0393_BHCLTKDSX2/FASTQ/"
INPUT_DIR_2=$INPUT_DIR_b1s2"210721_A00711_0400_BHF2YWDSX2/FASTQ/"
INPUT_DIR_3=$INPUT_DIR_b1s2"210728_A00711_0405_BHF5JLDSX2/FASTQ/"
INPUT_DIR_4=$INPUT_DIR_b1s2"210901_A00711_0420_AHGJJLDSX2/FASTQ/"

## results
OUTPUT_DIR=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/04_cellranger_count/02_set2"

## original sample IDs, correpond to FASTQ file names
SAMPLE_ID_FILE=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/90_sample_info/samples_set2_4_runs.txt"

## corresponding new sample names to be set by --id
SAMPLE_ID_NEW_FILE=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/90_sample_info/samples_set2_rename.txt"

## CR reference
## updated CR reference based on ensembl v108 (Gencode v42)
## from get_cellranger_reference.sh
CR_REF=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/03_cellranger_reference/GRCh38"


#-----------------------------------------------------------------------
# Cell Ranger count

mkdir -p $OUTPUT_DIR

N=${SLURM_ARRAY_TASK_ID}

SAMPLE_ID=$(cat $SAMPLE_ID_FILE | tail -n+${N} | head -1)
SAMPLE_ID_NEW=$(cat $SAMPLE_ID_NEW_FILE | tail -n+${N} | head -1)

echo "Input dir(s): "$INPUT_DIR_1", "$INPUT_DIR_2", "$INPUT_DIR_3", "$INPUT_DIR_4
echo "Sample: "$SAMPLE_ID
echo "Sample new: "$SAMPLE_ID_NEW

cd $OUTPUT_DIR

if [[ "$SAMPLE_ID_NEW" == "V_19" ]]
then
	$CELL_RANGER count --id=$SAMPLE_ID_NEW \
	--fastqs=$INPUT_DIR_1 \
	--sample=$SAMPLE_ID \
	--transcriptome=$CR_REF \
	--localcores=40 \
	--localmem=330 \
	--include-introns=true
elif [[ "$SAMPLE_ID_NEW" == "V_20" ]]
then
	$CELL_RANGER count --id=$SAMPLE_ID_NEW \
	--fastqs=$INPUT_DIR_1,$INPUT_DIR_3,$INPUT_DIR_4 \
	--sample=$SAMPLE_ID \
	--transcriptome=$CR_REF \
	--localcores=40 \
	--localmem=330 \
	--include-introns=true
else
	$CELL_RANGER count --id=$SAMPLE_ID_NEW \
	--fastqs=$INPUT_DIR_1,$INPUT_DIR_2,$INPUT_DIR_3,$INPUT_DIR_4 \
	--sample=$SAMPLE_ID \
	--transcriptome=$CR_REF \
	--localcores=40 \
	--localmem=330 \
	--include-introns=true
fi

## optional:
## R2 length is longer (151bp) than usually recommended/required (90bp)
## --r2-length=90

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
