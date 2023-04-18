#!/bin/bash

#SBATCH -p c_highmem_dri1
#SBATCH --job-name=FQ
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-144%144
#SBATCH --mem-per-cpu=5000 # memory limit per core
#####  #SBATCH --mem=96000 # memory limit per compute node for the job
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

## FASTQ files
## This is the directory for batch 2
INPUT_DIR="/gluster/dri02/rdsmbh/shared/rdsmbh/230327_A00748_0368_AH5CTMDSX5_fastq/"
## This if for batch 1 set 1
INPUT_DIR_b1s1="/gluster/dri02/rdscw/shared/webber/Endo_10X/FASTQ/Set_1/211008_A00748_0157_AHT7TJDSX2_fastq_L2_3_4/"
## This is for batch 1 set 2
INPUT_DIR_b1s2="/gluster/dri02/rdscw/shared/webber/Endo_10X/FASTQ/Set_2/"

OUTPUT_DIR="/scratch/c.mpmgb/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/01_fastqc/01_set1"
OUTPUT_DIR2="/scratch/c.mpmgb/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/01_fastqc/02_set2"
OUTPUT_DIR3="/scratch/c.mpmgb/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/01_fastqc/03_set3"

## Load FastQC module
module load FastQC

#-----------------------------------------------------------------------

mkdir -p $OUTPUT_DIR

N=${SLURM_ARRAY_TASK_ID}

FASTQ_FILE=$(find $INPUT_DIR -mindepth 1 -maxdepth 1 -name '*_001.fastq.gz' | sort | tail -n+${N} | head -1)
fastqc -o $OUTPUT_DIR -f fastq --noextract --quiet -t 1 $FASTQ_FILE

FASTQ_FILE=$(find $INPUT_DIR_b1s1 -mindepth 1 -maxdepth 1 -name '*_001.fastq.gz' | sort | tail -n+${N} | head -1)
fastqc -o $OUTPUT_DIR2 -f fastq --noextract --quiet -t 1 $FASTQ_FILE

FASTQ_FILE=$(find $INPUT_DIR_b1s2 -mindepth 1 -maxdepth 3 -name '*_001.fastq.gz' | sort | tail -n+${N} | head -1)
fastqc -o $OUTPUT_DIR3 -f fastq --noextract --quiet -t 1 $FASTQ_FILE

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
