#!/bin/bash
#SBATCH --job-name=MAGMA_step3
#SBATCH -p c_compute_dri1
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --time=0-01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --account=scw1329
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bernardo-harringtong@cardiff.ac.uk

#####SBATCH --mem=120G

echo "*****************************************************************"
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
module load parallel


WORK_DIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/994_magma_inputs/"
cd $WORK_DIR
export WORK_DIR

# Populate Set_Annot_Files array with files matching '*.magma.txt'
Set_Annot_Files=(${WORK_DIR}/*.magma.txt)
# Populate Gene_Results_Files array with files matching '*.genes.raw'
Gene_Results_Files=(${WORK_DIR}/gwas_backgrounds/*.genes.raw)

# Function to run magma with given files and output prefix
run_magma() {
    local Set_Annot_File=$1
    local Gene_Results_File=$2
    local annot_name=$(basename "${Set_Annot_File}" .magma.txt)
    local gene_name=$(basename "${Gene_Results_File}" .genes.raw)
    local Output_Prefix="./results/${annot_name}_${gene_name}"

    # Run magma command
    magma \
        --gene-results $Gene_Results_File \
        --set-annot $Set_Annot_File \
        --out $Output_Prefix
}

# Export the function so it can be used by parallel
export -f run_magma

# Use parallel to run the combinations
parallel run_magma ::: "${Set_Annot_Files[@]}" ::: "${Gene_Results_Files[@]}"


echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"