#!/bin/bash

#SBATCH -p c_vhighmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=bbb_level1_ldsc_25k
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-22
##### #SBATCH --mem-per-cpu=8000 # memory limit per core
#SBATCH --mem=300G # memory limit per compute node for the job
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


N=${SLURM_ARRAY_TASK_ID}


LDSC="/LDSC/ldsc"
#GENESET="/TG3_level1_GeneSet"
PHASE3="/nfs/dri/02/rdscw/shared/public/LDSc_regression
/1000G_EUR_Phase3_plink"
OUTPUT="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/995_ldsc_inputs/03_annotations"
ENSG ="/nfs/dri/02/rdscw/shared/public/LDSc_regression/"
INPUT_FILES="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/995_ldsc_inputs/"

mkdir -p $OUTPUT

# for file in Emeka_Fibroblast Emeka_Immune Emeka_Satglia Emeka_Schwann_M Emeka_Vascular Emeka_control_Fibroblast Emeka_control_Immune Emeka_control_Satglia Emeka_control_Schwann_M Emeka_control_Vascular Emeka_control_neuronal Emeka_neuronal
# do
# echo ${file}
# 
#   python $LDSC/make_annot.py \
#   --gene-set-file $GENESET/$file.GeneSet \
#   --gene-coord-file $ENSG/ENSG_coord.txt \
#   --windowsize 25000 \
#   --bimfile $PHASE3/1000G.EUR.QC.$N.bim \
#   --annot-file $OUTPUT/$file.$N.annot.gz
# done

for file in ${INPUT_FILE}*; do
  basefile=$(basename "$file")
  echo "${basefile}"

  python $LDSC/make_annot.py \
  --gene-set-file $file \
  --gene-coord-file $ENSG/ENSG_coord.txt \
  --windowsize 25000 \
  --bimfile $PHASE3/1000G.EUR.QC.$N.bim \
  --annot-file $OUTPUT/${basefile}.$N.annot.gz
done


echo -e "\n************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
