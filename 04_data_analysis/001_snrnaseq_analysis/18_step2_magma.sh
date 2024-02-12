#!/bin/bash
#SBATCH --job-name=MAGMA_step2
#SBATCH -p c_highmem_dri1
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --time=0-10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --array=1
#SBATCH --mem=120G
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
module load parallel

WORK_DIR="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/993_magma_inputs/"
cd $WORK_DIR

#### Without proxies and with APOE - 35k10k window

mkdir temp_annot_35k10k # make a temporary directory to host the intermediate files

# PLINK files from reference (.bed/.bim/.fam)
Data_File=${WORK_DIR}"g1000_eur"
# Output of step 1
Annot_File=${WORK_DIR}"NCBI37_annotated_window_35k10k.genes.annot"
# GWAS summary stats
# Ata has summary stats here: /scratch/scw1751/shared/Bellenguez_2022_sumstats/
# GCST90027158_buildGRCh38_noATGC.tsv is for build 38, but the MAGMA site has reference data for build 37 only, asking sam for link to 38
# Should try with and without APOE in summary stats
# Should also try with and without proxy cases, file: EADB_release_Feb2020.meta_model_pcs.mac_info_20_GC_OFF.formatted.exc_perc_cases_50.het_5e-8.freq_amp_40_noATGC_GRCh37.tsv
SNP_Pval_File="/scratch/scw1751/shared/Bellenguez_2022_sumstats/EADB_release_Feb2020.meta_model_pcs.mac_info_20_GC_OFF.formatted.exc_perc_cases_50.het_5e-8.freq_amp_40_noATGC_GRCh37.tsv"

Output_Prefix="EUROPE_noproxy_35k10k"

    # run magma in parallel, 10 threads in this case

parallel magma \
   --batch {} 10 \
   --bfile $Data_File \
   --gene-annot $Annot_File \
   --gene-model snp-wise=mean \
   --pval $SNP_Pval_File ncol=N \
   --out temp_annot_35k10k/$Output_Prefix \
::: {1..10}

# merge all intermediate files generated under the temp_annot files
# and send out for one single file set

magma \
   --merge temp_annot_35k10k/$Output_Prefix \
   --out temp_annot_35k10k/$Output_Prefix

# extract merged files for subsequent analysis

cp ./temp_annot_35k10k/$Output_Prefix.genes.* .

# remove the temporary directory

rm -r temp_annot_35k10k

#### Without proxies and with APOE - no window

mkdir temp_annot # make a temporary directory to host the intermediate files

Annot_File=${WORK_DIR}"NCBI37_annotated_nowindow.genes.annot"
SNP_Pval_File="/scratch/scw1751/shared/Bellenguez_2022_sumstats/EADB_release_Feb2020.meta_model_pcs.mac_info_20_GC_OFF.formatted.exc_perc_cases_50.het_5e-8.freq_amp_40_noATGC_GRCh37.tsv"
Output_Prefix="EUROPE_noproxy_nowindow"

# run magma in parallel, 10 threads in this case
parallel magma \
   --batch {} 10 \
   --bfile $Data_File \
   --gene-annot $Annot_File \
   --gene-model snp-wise=mean \
   --pval $SNP_Pval_File ncol=N \
   --out temp_annot/$Output_Prefix \
::: {1..10}

# merge all intermediate files generated under the temp_annot files
# and send out for one single file set

magma \
   --merge temp_annot/$Output_Prefix \
   --out temp_annot/$Output_Prefix

# extract merged files for subsequent analysis

cp ./temp_annot/$Output_Prefix.genes.* .

# remove the temporary directory

rm -r temp_annot


#### With proxies and APOE - 35k10k window

mkdir temp_annot_35k10k_nobb # make a temporary directory to host the intermediate files

Annot_File=${WORK_DIR}"NCBI37_annotated_window_35k10k.genes.annot"
SNP_Pval_File="/scratch/scw1751/shared/Bellenguez_2022_sumstats/GCST90027158_buildGRCh37_noATGC_forMAGMA.tsv"
Output_Prefix="EUROPEUKBB_35k10k"


# run magma in parallel, 10 threads in this case
parallel magma \
   --batch {} 10 \
   --bfile $Data_File \
   --gene-annot $Annot_File \
   --gene-model snp-wise=mean \
   --pval $SNP_Pval_File ncol=N \
   --out temp_annot_35k10k_nobb/$Output_Prefix \
::: {1..10}

# merge all intermediate files generated under the temp_annot files
# and send out for one single file set

magma \
   --merge temp_annot_35k10k_nobb/$Output_Prefix \
   --out temp_annot_35k10k_nobb/$Output_Prefix

# extract merged files for subsequent analysis

cp ./temp_annot_35k10k_nobb/$Output_Prefix.genes.* .

# remove the temporary directory

rm -r temp_annot_35k10k_nobb

#### With proxies and APOE - no window

mkdir temp_annot_nobb # make a temporary directory to host the intermediate files

Annot_File=${WORK_DIR}"NCBI37_annotated_nowindow.genes.annot"
SNP_Pval_File="/scratch/scw1751/shared/Bellenguez_2022_sumstats/GCST90027158_buildGRCh37_noATGC_forMAGMA.tsv"
Output_Prefix="EUROPEUKBB_nowindow"


# run magma in parallel, 10 threads in this case

parallel magma \
   --batch {} 10 \
   --bfile $Data_File \
   --gene-annot $Annot_File \
   --gene-model snp-wise=mean \
   --pval $SNP_Pval_File ncol=N \
   --out temp_annot_nobb/$Output_Prefix \
::: {1..10}

# merge all intermediate files generated under the temp_annot files
# and send out for one single file set

magma \
   --merge temp_annot_nobb/$Output_Prefix \
   --out temp_annot_nobb/$Output_Prefix

# extract merged files for subsequent analysis

cp ./temp_annot_nobb/$Output_Prefix.genes.* .

# remove the temporary directory

rm -r temp_annot_nobb

#### With proxies and without APOE - no window

mkdir temp_annot_nobb # make a temporary directory to host the intermediate files

Annot_File=${WORK_DIR}"NCBI37_annotated_nowindow.genes.annot"
SNP_Pval_File="/scratch/scw1751/shared/Bellenguez_2022_sumstats/GCST90027158_buildGRCh37_noATGC_noAPOE_forMAGMA.tsv"
Output_Prefix="EUROPEUKBB_nowindow_noAPOE"

# run magma in parallel, 10 threads in this case
parallel magma \
   --batch {} 10 \
   --bfile $Data_File \
   --gene-annot $Annot_File \
   --gene-model snp-wise=mean \
   --pval $SNP_Pval_File ncol=N \
   --out temp_annot_nobb/$Output_Prefix \
::: {1..10}

# merge all intermediate files generated under the temp_annot files
# and send out for one single file set

magma \
   --merge temp_annot_nobb/$Output_Prefix \
   --out temp_annot_nobb/$Output_Prefix

# extract merged files for subsequent analysis

cp ./temp_annot_nobb/$Output_Prefix.genes.* .

# remove the temporary directory

rm -r temp_annot_nobb

#### With proxies and without APOE - 35k10k window

mkdir temp_annot_nobb # make a temporary directory to host the intermediate files

Annot_File=${WORK_DIR}"NCBI37_annotated_window_35k10k.genes.annot"
SNP_Pval_File="/scratch/scw1751/shared/Bellenguez_2022_sumstats/GCST90027158_buildGRCh37_noATGC_noAPOE_forMAGMA.tsv"
Output_Prefix="EUROPEUKBB_35k10k_noAPOE"

# run magma in parallel, 10 threads in this case
parallel magma \
   --batch {} 10 \
   --bfile $Data_File \
   --gene-annot $Annot_File \
   --gene-model snp-wise=mean \
   --pval $SNP_Pval_File ncol=N \
   --out temp_annot_nobb/$Output_Prefix \
::: {1..10}

# merge all intermediate files generated under the temp_annot files
# and send out for one single file set

magma \
   --merge temp_annot_nobb/$Output_Prefix \
   --out temp_annot_nobb/$Output_Prefix

# extract merged files for subsequent analysis

cp ./temp_annot_nobb/$Output_Prefix.genes.* .

# remove the temporary directory

rm -r temp_annot_nobb

echo -e "\n*****************************************************************"
echo "Finished at: "`date`
echo "*****************************************************************"