#!/bin/bash

#SBATCH -p c_highmem_dri1
#SBATCH --job-name=scflow-seurat-processing
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G # memory limit per compute node for the job
#SBATCH --time=1-00:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --mail-user Bernardo-HarringtonG@cardiff.ac.uk # email on fail
#SBATCH --mail-type END,FAIL

# load environment
# note mamba env created with the following:
#mamba create --prefix r-env -c conda-forge -c bioconda -c R r-base=4.2 r-seurat r-future

# activate conda environment
module load anaconda
source activate

# move to directory
cd /scratch/scw1329/gmbh/blood-brain-barrier-in-ad

conda activate /scratch/scw1329/gmbh/blood-brain-barrier-in-ad/r-env

Rscript 04_data_analysis/990_code_libraries/scflow-seurat-processing.R
