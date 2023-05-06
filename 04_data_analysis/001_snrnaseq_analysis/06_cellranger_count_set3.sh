#!/bin/bash

#SBATCH -p c_vhighmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=cellranger_count_set3
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-16%2
##### #SBATCH --mem-per-cpu=8000 # memory limit per core
#SBATCH --mem=740G # memory limit per compute node for the job
#SBATCH --time=4-00:00 # maximum job time in D-HH:MM
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


## FASTQ files, set 2 with 16 samples, data from 4 runs, sequenced in Oxford
## NOTE, V_19, only sequenced in 1 run (INPUT_DIR_1)
## NOTE, V_20, 1 run (INPUT_DIR_2) only with a single read, skip that run
## handled via if statement below
INPUT_DIR_1="/gluster/dri02/rdsmbh/shared/rdsmbh/230327_A00748_0368_AH5CTMDSX5_fastq/"

## results
OUTPUT_DIR="/scratch/c.mpmgb/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/04_cellranger_count/03_set3"

## original sample IDs, correpond to FASTQ file names
SAMPLE_ID_FILE="/scratch/c.mpmgb/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/90_sample_info/samples_set3.txt"

## CR executable
CELL_RANGER="/scratch/c.mpmgb/tools/cellranger-7.1.0/bin/cellranger"

## CR reference
## pre-computed and downloaded from CR website as is (refdata-gex-GRCh38-2020-A.tar.gz)
## based on ensembl v98/gencode v32
## CR_REF="/scratch/c.mpmfw/Tools/CellRanger/CellRanger_references/refdata-gex-GRCh38-2020-A"
## updated CR reference based on ensembl v108 (Gencode v42)
## from get_cellranger_reference.sh
CR_REF="/scratch/c.mpmgb/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/03_cellranger_reference/GRCh38"


#-----------------------------------------------------------------------
# Cell Ranger count

mkdir -p $OUTPUT_DIR

N=${SLURM_ARRAY_TASK_ID}

SAMPLE_ID=$(cat $SAMPLE_ID_FILE | tail -n+${N} | head -1)
# SAMPLE_ID_NEW=$(cat $SAMPLE_ID_NEW_FILE | tail -n+${N} | head -1)

echo "Input dir: "$INPUT_DIR
echo "Sample: "$SAMPLE_ID
# echo "Sample new: "$SAMPLE_ID_NEW

cd $OUTPUT_DIR

$CELL_RANGER count --id=$SAMPLE_ID \
--fastqs=$INPUT_DIR \
--sample=$SAMPLE_ID \
--transcriptome=$CR_REF \
--localcores=40 \
--localmem=600 \
--include-introns=true \
--no-bam

## optional:
## R2 length is longer (151bp) than usually recommended/required (90bp)
## --r2-length=90

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
