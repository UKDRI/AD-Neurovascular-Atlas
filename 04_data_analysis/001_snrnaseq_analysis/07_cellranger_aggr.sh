#!/bin/bash

#SBATCH -p c_vhighmem_dri1
#SBATCH --job-name=cellranger_aggr
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem=540G # memory limit per compute node for the job
#SBATCH --time=1-00:00 # maximum job time in D-HH:MM
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
INPUT_CSV="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/05_cellranger_aggr/cellranger_ag.csv"

## name for directory to save output
DIR_NAME="set_aggr"

## results
OUTPUT_DIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/05_cellranger_aggr"

## CR executable
CELL_RANGER="/scratch/c.mpmgb/tools/cellranger-7.1.0/bin/cellranger"


#-----------------------------------------------------------------------
# Cell Ranger aggr

mkdir -p $OUTPUT_DIR

cd $OUTPUT_DIR

$CELL_RANGER aggr --id=$DIR_NAME \
--csv=$INPUT_CSV

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
