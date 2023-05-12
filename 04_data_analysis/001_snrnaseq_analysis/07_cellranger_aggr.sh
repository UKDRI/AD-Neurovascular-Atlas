#!/bin/bash

#SBATCH -p c_vhighmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=cellranger_aggr
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-4%2
##### #SBATCH --mem-per-cpu=8000 # memory limit per core
#SBATCH --mem=640G # memory limit per compute node for the job
#SBATCH --time=2-10:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --mail-user Bernardo-HarringtonG@cardiff.ac.uk # email on fail
#SBATCH --mail-type END,FAIL

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


## path to input csv with samples IDs and paths to cellranger count outputs
INPUT_CSVS="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/05_cellranger_aggr/cellranger_aggr_csvs.txt"


## results
OUTPUT_DIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/05_cellranger_aggr"

## CR executable
CELL_RANGER="/scratch/c.mpmgb/tools/cellranger-7.1.0/bin/cellranger"


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
