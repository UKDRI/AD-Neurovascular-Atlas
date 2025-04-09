#!/bin/bash

#SBATCH -p c_vhighmem_dri1
#SBATCH --job-name=run_tradeseq_split
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=720G # memory limit per compute node for the job
#SBATCH --time=3-00:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%A_%a_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%A_%a_%J.txt
#SBATCH --mail-user Bernardo-HarringtonG@cardiff.ac.uk # email on fail
#SBATCH --mail-type END,FAIL

# load environment
# note mamba env created with the following:
#mamba create --prefix r-env -c conda-forge -c bioconda -c R -c paul.martin-2 r-base r-future r-readr r-remotes r-here r-tidyverse r-qs bioconductor-slingshot bioconductor-tradeseq bioconductor-trajectoryutils bioconductor-delayedmatrixstats

# activate conda environment
module load anaconda
source activate

# move to directory
cd /scratch/scw1329/gmbh/blood-brain-barrier-in-ad

conda activate /scratch/scw1329/gmbh/blood-brain-barrier-in-ad/r-env

Rscript 04_data_analysis/990_code_libraries/run_tradeseq_split.R
