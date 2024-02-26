
#!/bin/bash

#SBATCH -p c_vhighmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=TG3_level1_ldsc_25k
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-22
##### #SBATCH --mem-per-cpu=8000 # memory limit per core
#SBATCH --mem=300G # memory limit per compute node for the job
#SBATCH --time=4-00:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/.../%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/.../%x_err_%A_%a_%J.txt

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
GWASSTAT="/TG3_level1_GeneSet"
BASELINE="/1000G_EUR_Phase3_baseline"
OUTPUT="path"
LDCT="/path"
WEIGHTS="/weights_hm3_no_hla"


for gwas in multi_sited neuro_paindoecho ${gwas}python2 $LDSC/ldsc.py \--h2-cts $GWASSTAT/$gwas.sumstats.gz \--ref-ld-chr $BASELINE/baseline. \--out $OUTPUT/$gwas.10.25k \--ref-ld-chr-cts $LDCT/TG10_level1_25k_1000Gv3.ldcts \--w-ld-chr $WEIGHTS/weights. done

echo -e "\n************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"
