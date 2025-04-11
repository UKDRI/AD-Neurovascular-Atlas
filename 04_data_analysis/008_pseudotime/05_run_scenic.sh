#!/bin/bash
#SBATCH -p c_highmem_dri1
#SBATCH --job-name=run_scenic_endomt
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=350G # memory limit per compute node for the job
#SBATCH --time=3-00:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --mail-user Bernardo-HarringtonG@cardiff.ac.uk # email on fail
#SBATCH --mail-type END,FAIL

# load environment
# note mamba env created with the following:
#cd /scratch/scw1329/gmbh/blood-brain-barrier-in-ad
#mamba create --prefix pyscenic pyscenic python=3.10

# load environment
module load anaconda
source activate

# move to directory
cd /scratch/scw1329/gmbh/blood-brain-barrier-in-ad
conda activate /scratch/scw1329/gmbh/blood-brain-barrier-in-ad/pyscenic

DB_DIR=/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/990_processed_data/008_pseudotime

# Step 1: Run the Python script for GRNBoost2 network inference
echo "Starting GRNBoost2 network inference..."
python run_grnboost2.py \
    --input ${DB_DIR}/scenic_input_counts.csv \
    --tf-list ${DB_DIR}/hs_hgnc_tfs.txt \
    --output ${DB_DIR}/grnboost2_network.tsv \
    --cores 30

# Step 2: Run pyscenic ctx for cisTarget analysis
echo "Starting cisTarget analysis..."
pyscenic ctx ${DB_DIR}/grnboost2_network.tsv \
  ${DB_DIR}/01_scenic/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather,\
${DB_DIR}/01_scenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
  --annotations_fname ${DB_DIR}/01_scenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname ${DB_DIR}/scenic_input_counts.csv \
  --output ${DB_DIR}/reg.csv \
  --mode "dask_multiprocessing" \
  --num_workers 30

# Step 3: Run pyscenic aucell for cellular enrichment
echo "Starting AUCell analysis..."
pyscenic aucell ${DB_DIR}/scenic_input_counts.csv \
    ${DB_DIR}/reg.csv \
    --output ${DB_DIR}/auc_mtx.csv \
    --num_workers 30

echo "SCENIC analysis completed."
