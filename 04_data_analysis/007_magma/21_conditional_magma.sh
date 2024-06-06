#!/bin/bash

module purge
module load magma/1.10
module load parallel


WORK_DIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/994_magma_inputs/"
cd $WORK_DIR
export WORK_DIR

# Populate Set_Annot_Files array with files matching '*.magma.txt'
Set_Annot_Files=magma_celltypes_level2_all_controls_top-ten-percent.magma.txt
Set_Annot_File=magma_celltypes_level2_all_controls_top-ten-percent.magma.txt
export Set_Annot_File

# Populate Gene_Results_Files array with files matching '*.genes.raw'
Gene_Results_Files=gwas_backgrounds/EUROPEUKBB_35k10k.genes.raw
Gene_Results_File=gwas_backgrounds/EUROPEUKBB_35k10k.genes.raw
export Gene_Results_File

# Add the conditional gene set files
Conditional_Geneset_Files=(${WORK_DIR}/sig_genelist/*.condition.txt)
conditional_array=("Microglia-activated" "Pericyte-2" "Perivascular-FB-KAZN2")

Output_Prefix="./sig_genelist/results/level2_all_controls_top-ten-percent_AD_microglia-activated_conditional"

    magma \
        --gene-results $Gene_Results_Files \
        --set-annot $Set_Annot_Files \
        --model condition=${conditional_array[0]}\
        --out $Output_Prefix
        
Output_Prefix="./sig_genelist/results/level2_all_controls_top-ten-percent_AD_pericyte-2_conditional"
    magma \
        --gene-results $Gene_Results_Files \
        --set-annot $Set_Annot_Files \
        --model condition=${conditional_array[1]}\
        --out $Output_Prefix
        
Output_Prefix="./sig_genelist/results/level2_all_controls_top-ten-percent_AD_perivascular-fb-kazn2_conditional"
    magma \
        --gene-results $Gene_Results_Files \
        --set-annot $Set_Annot_Files \
        --model condition=${conditional_array[2]}\
        --out $Output_Prefix
        
# Function to run magma with given files and output prefix
run_magma() {
    local Conditional=$1
    local annot_name=$(basename "${Set_Annot_File}" .magma.txt)
    local gene_name=$(basename "${Gene_Results_File}" .genes.raw)
    local Output_Prefix="./results/${annot_name}_${gene_name}_${Conditional}"

    # Run magma command
    magma \
        --gene-results $Gene_Results_File \
        --set-annot $Set_Annot_File \
        --model condition=$Conditional \
        --out $Output_Prefix
}

# Export the function so it can be used by parallel
export -f run_magma

# Use parallel to run the combinations
parallel run_magma ::: "${conditional_array[@]}"
