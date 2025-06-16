#!/bin/bash

#SBATCH -p c_highmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=cellranger_count_set3
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-40%14
#SBATCH --mem=340G # memory limit per compute node for the job
#SBATCH --time=2-10:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt

## Project root variables defined in the following
source /scratch/c.mpmgb/blood-brain-barrier-in-ad/04_data_analysis/001_snrnaseq_analysis/00_hpc_variables.sh

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

## results
OUTPUT_DIR=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/04_cellranger_count/03_set3"

## original sample IDs, correpond to FASTQ file names
SAMPLE_ID_FILE=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/90_sample_info/samples_set3.txt"

## CR reference
## updated CR reference based on ensembl v108 (Gencode v42)
## from get_cellranger_reference.sh
CR_REF=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/03_cellranger_reference/GRCh38"

#-----------------------------------------------------------------------
# Cell Ranger count

mkdir -p $OUTPUT_DIR

N=${SLURM_ARRAY_TASK_ID}

SAMPLE_ID=$(cat $SAMPLE_ID_FILE | tail -n+${N} | head -1)
# Get the sample ID string after the first "_"
SAMPLE_ID_NEW="${SAMPLE_ID#*_}"

echo "Input dir: "$INPUT_DIR
echo "Sample: "$SAMPLE_ID
echo "Sample new: "$SAMPLE_ID_NEW

cd $OUTPUT_DIR

$CELL_RANGER count --id=$SAMPLE_ID_NEW \
--fastqs=$INPUT_DIR \
--sample=$SAMPLE_ID \
--transcriptome=$CR_REF \
--localcores=40 \
--localmem=330 \
--include-introns=true

## optional:
## R2 length is longer (151bp) than usually recommended/required (90bp)
## --r2-length=90

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
