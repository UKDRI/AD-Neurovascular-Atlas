#!/bin/bash

#SBATCH -p c_compute_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=bbb_25K_ldscore_regression
#SBATCH --mem=100G # memory limit per compute node for the job
#SBATCH --time=0-10:00 # maximum job time in D-HH:MM
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
module load parallel

## NOTE: I think maybe the .annot files need to be in the same dir
# code to remove .tsv from file names in case it matters
# for file in *tsv*
# do
#     # Remove the .tsv part from the filename
#     new_name="${file/.tsv/}"
# 
#     # Optionally, you can rename the files if that's the intention
#     mv "$file" "$new_name"
# 
#     echo "Original filename: $file"
#     echo "New filename: $new_name"
# done

# # Loop over all files that contain a hyphen
# for file in *-*; do
#     # Replace all hyphens with underscores in the filename
#     new_name="${file//-/_}"
# 
#     # Rename the file
#     mv "$file" "$new_name"
# 
#     # Echo the old and new filenames (optional)
#     echo "Renamed $file to $new_name"
# done

BASEDIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/995_ldsc_inputs"
LDSC=$BASEDIR"/ldsc"
BASELINE=$BASEDIR"/1000G_EUR_Phase3_baseline"
WEIGHTS=$BASEDIR"/weights_hm3_no_hla"
OUTPUT=$BASEDIR"/05_ldscore_regression"
LDCT=$BASEDIR"/celltype_markers.ldcts"

mkdir -p $OUTPUT

GWASSTAT=$BASEDIR"/gwas_backgrounds"

# Loop over all .sumstats.gz files in the directory
for gwas_file in "$GWASSTAT"/*.sumstats.gz
do
    # Extract the base name without the extension for use in the output file name
    gwas=$(basename "$gwas_file" .sumstats.gz)

    echo "${gwas}"

    python2 $LDSC/ldsc.py \
    --h2-cts "$gwas_file" \
    --ref-ld-chr $BASELINE/baseline. \
    --out $OUTPUT/"$gwas".25k \
    --ref-ld-chr-cts $LDCT \
    --w-ld-chr $WEIGHTS/weights.
done

echo -e "\n************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"