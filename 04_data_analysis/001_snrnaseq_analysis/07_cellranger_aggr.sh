#!/bin/bash

#SBATCH -p c_vhighmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=cellranger_aggr
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-4%2
#SBATCH --mem=640G # memory limit per compute node for the job
#SBATCH --time=2-10:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt

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


## path to input csvs with samples IDs and paths to cellranger count outputs
INPUT_CSVS=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/05_cellranger_aggr/cellranger_aggr_csvs.txt"

## results
OUTPUT_DIR=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/05_cellranger_aggr"

#-----------------------------------------------------------------------
# Cell Ranger aggr

mkdir -p $OUTPUT_DIR

N=${SLURM_ARRAY_TASK_ID}

INPUT_CSV=$(cat $INPUT_CSVS | tail -n+${N} | head -1)
# Subset to last part of the file path
INPUT_CSV_NAME="${INPUT_CSV##*/}"
# Subset string to first "."
SAVE_DIR="${INPUT_CSV_NAME%%.*}"

echo "Input CSV: "$INPUT_CSV
echo "Save dir: "$SAVE_DIR

cd $OUTPUT_DIR

$CELL_RANGER aggr --id=$SAVE_DIR \
--csv=$INPUT_CSV

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
