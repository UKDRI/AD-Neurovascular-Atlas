#!/bin/bash

#SBATCH -p c_vhighmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=bbb_25K_ldscore
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=30
#SBATCH --array=1-22
##### #SBATCH --array=1-7590%14
##### #SBATCH --mem-per-cpu=8000 # memory limit per core
#SBATCH --mem=360G # memory limit per compute node for the job
#SBATCH --time=3-00:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bernardo-harringtong@cardiff.ac.uk

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

module load ldsc/1.0.1
module load bedtools/2.29.2

## input
## gene sets

chrom=${SLURM_ARRAY_TASK_ID}
echo "Processing chromosome ${chrom}"

#echo "Input dir: "$INPUT_DIR

LDSC="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/995_ldsc_inputs/ldsc"
PHASE3="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/995_ldsc_inputs/1000G_EUR_Phase3_plink"
ANNOT="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/995_ldsc_inputs/03_annotations"
OUTPUT="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/995_ldsc_inputs/04_ldscores"
HAPMAP3="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/995_ldsc_inputs/hapmap3_snps"

for filepath in ${ANNOT}/*.tsv*; do  # Assuming the files have .tsv followed by more text
    file_with_path="${filepath%.tsv*}"  # Removes everything after .tsv, including .tsv itself
    basefile=$(basename -- "$file_with_path")  # Extracts the base filename without the directory path

    echo ${basefile}
  
  python $LDSC/ldsc.py \
  --l2 \
  --bfile $PHASE3/1000G.EUR.QC.$chrom \
  --ld-wind-cm 1 \
  --annot ${file_with_path}.tsv.$chrom.annot.gz \
  --thin-annot \
  --out $OUTPUT/${basefile}.$chrom \
  --print-snps $HAPMAP3/hm.$chrom.snp
done


echo -e "\n************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
