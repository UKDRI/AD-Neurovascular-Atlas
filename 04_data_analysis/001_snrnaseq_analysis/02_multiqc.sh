#!/bin/bash

#SBATCH -p c_highmem_dri1
#SBATCH --job-name=multiqc_bbb
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00 # maximum job time in D-HH:MM
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

# Load module
module load multiqc

## collect MultiQC reports
OUTPUT_DIR=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/02_multiqc/"

#-----------------------------------------------------------------------
# FastQC runs

## Endo 10X Vascular and Parenchymal fraction set 1 - 24 samples
INPUT_DIR=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/01_fastqc/01_set1/"
## Endo 10X Vascular and Parenchymal fraction set 2 - 16 samples
INPUT_DIR2=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/01_fastqc/02_set2/"
## Endo 10X Vascular and Parenchymal fraction set 3 - 40 samples
INPUT_DIR3=$PROJECT_ROOT"03_data/990_processed_data/001_snrnaseq/01_fastqc/03_set3/"

multiqc -o $OUTPUT_DIR"01_set1_reads" \
--filename "fastqc_endo_10X_set1_reads" \
--title "FastQC endo 10X set1 reads" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR""*_R[1-2]_*_fastqc.zip

multiqc -o $OUTPUT_DIR"01_set1_index" \
--filename "fastqc_endo_10X_set1_index" \
--title "FastQC endo 10X set1 index" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR""*_I[1-2]_*_fastqc.zip

multiqc -o $OUTPUT_DIR"02_set2_reads" \
--filename "fastqc_endo_10X_set2_reads" \
--title "FastQC endo 10X set2 reads" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR2""*_R[1-2]_*_fastqc.zip

multiqc -o $OUTPUT_DIR"02_set2_index" \
--filename "fastqc_endo_10X_set2_index" \
--title "FastQC endo 10X set2 index" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR2""*_I[1-2]_*_fastqc.zip

multiqc -o $OUTPUT_DIR"03_set3_reads" \
--filename "fastqc_endo_10X_set3_reads" \
--title "FastQC endo 10X set3 reads" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR3""*_R[1-2]_*_fastqc.zip

multiqc -o $OUTPUT_DIR"03_set3_index" \
--filename "fastqc_endo_10X_set3_index" \
--title "FastQC endo 10X set3 index" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR3""*_I[1-2]_*_fastqc.zip

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
