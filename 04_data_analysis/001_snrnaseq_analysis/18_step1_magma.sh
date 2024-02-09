#!/bin/bash
#SBATCH --job-name=MAGMA_step1
#SBATCH -p c_highmem_dri1
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --time=0-05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1
#SBATCH --account=scw1329
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bernardo-harringtong@cardiff.ac.uk

echo "*****************************************************************"
echo "All jobs in this array have:"
echo "  - SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID}"
echo "  - SLURM_ARRAY_TASK_COUNT: ${SLURM_ARRAY_TASK_COUNT}"
echo "  - SLURM_ARRAY_TASK_MIN: ${SLURM_ARRAY_TASK_MIN}"
echo "  - SLURM_ARRAY_TASK_MAX: ${SLURM_ARRAY_TASK_MAX}"
echo "Job in the array has:"
echo "    - SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "    - SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "Host: "`hostname`
echo "Number of threads (nproc): "`nproc`
echo "Total memory in GB: "`free -g | grep -oP '\d+' | sed -n 1p`
echo "Used memory in GB: "`free -g | grep -oP '\d+' | sed -n 2p`
echo "Free memory in GB: "`free -g | grep -oP '\d+' | sed -n 3p`
echo "Username: "`whoami`
echo "Started at: "`date`
echo -e "*****************************************************************\n"

module purge
module load magma/1.10

WORK_DIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/993_magma_inputs/"
cd $WORK_DIR

# SNP file from the reference
SNP_Loc_File=${WORK_DIR}"g1000_eur.bim"
# Gene location file - make sure build is correct
Gene_Loc_File=${WORK_DIR}"NCBI37.3.gene.loc"
Output_Prefix="NCBI37_annotated_window"

magma \
    --annotate window=35,10 \
    --snp-loc $SNP_Loc_File \
    --gene-loc $Gene_Loc_File \
    --out $Output_Prefix

Output_Prefix="NCBI37_annotated_nowindow"

magma \
    --annotate window=0\
    --snp-loc $SNP_Loc_File \
    --gene-loc $Gene_Loc_File \
    --out $Output_Prefix


echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"