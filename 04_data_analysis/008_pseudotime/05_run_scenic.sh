#!/bin/bash
#SBATCH -p c_highmem_dri1
#SBATCH --job-name=run_scenic_endomt
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=350G
#SBATCH --time=3-00:00
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --mail-user=Bernardo-HarringtonG@cardiff.ac.uk
#SBATCH --mail-type=END,FAIL

# Load environment
# note mamba env created with the following:
#cd /scratch/scw1329/gmbh/blood-brain-barrier-in-ad
#mamba create --prefix pyscenic pyscenic python=3.10
module load anaconda
source activate
cd /scratch/scw1329/gmbh/blood-brain-barrier-in-ad
conda activate /scratch/scw1329/gmbh/blood-brain-barrier-in-ad/pyscenic

DB_DIR=/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/990_processed_data/008_pseudotime
OUT_DIR=${DB_DIR}

# Check input files exist
COUNTS_FILE=${DB_DIR}/scenic_input_counts.csv
TF_LIST=${DB_DIR}/hs_hgnc_tfs.txt
MOTIF_ANNOT=${DB_DIR}/01_scenic/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
DATABASE1=${DB_DIR}/01_scenic/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
DATABASE2=${DB_DIR}/01_scenic/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

# Function for checking file existence
check_file () {
  if [ ! -f "$1" ]; then
    echo "ERROR: Required file $1 not found. Exiting."
    exit 1
  fi
}

# Check required input files
echo "Checking required input files..."
check_file "$COUNTS_FILE"
check_file "$TF_LIST"
check_file "$MOTIF_ANNOT"
check_file "$DATABASE1"
check_file "$DATABASE2"

# Step 1: Run GRNBoost2 if output doesn't exist
GRNBOOST_OUTPUT=${OUT_DIR}/grnboost2_network.tsv
if [ -f "$GRNBOOST_OUTPUT" ]; then
  echo "GRNBoost2 output already exists at $GRNBOOST_OUTPUT. Skipping GRNBoost2 step."
else
  echo "Starting GRNBoost2 network inference..."
  python /scratch/scw1329/gmbh/blood-brain-barrier-in-ad/04_data_analysis/008_pseudotime/04_scenic_python.py \
    --input "$COUNTS_FILE" \
    --tf-list "$TF_LIST" \
    --output "$GRNBOOST_OUTPUT" \
    --cores 30
fi

# Step 2: Run pyscenic ctx if output doesn't exist
REG_OUTPUT=${OUT_DIR}/reg.csv
if [ -f "$REG_OUTPUT" ]; then
  echo "cisTarget results already exist at $REG_OUTPUT. Skipping ctx step."
else
  echo "Starting cisTarget analysis..."
  pyscenic ctx "$GRNBOOST_OUTPUT" \
    "$DATABASE1" "$DATABASE2" \
    --annotations_fname "$MOTIF_ANNOT" \
    --expression_mtx_fname "$COUNTS_FILE" \
    --output "$REG_OUTPUT" \
    --mode "dask_multiprocessing" \
    --num_workers 30
fi

# Step 3: Run pyscenic aucell if output doesn't exist
AUCELL_OUTPUT=${OUT_DIR}/auc_mtx.csv
if [ -f "$AUCELL_OUTPUT" ]; then
  echo "AUCell output already exists at $AUCELL_OUTPUT. Skipping AUCell step."
else
  echo "Starting AUCell analysis..."
  pyscenic aucell "$COUNTS_FILE" \
    "$REG_OUTPUT" \
    --output "$AUCELL_OUTPUT" \
    --num_workers 30
fi

echo "SCENIC analysis completed."
