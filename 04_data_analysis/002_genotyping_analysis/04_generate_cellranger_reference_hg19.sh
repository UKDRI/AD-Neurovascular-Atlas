#!/bin/bash

#SBATCH -p c_highmem_dri1 	###c_highmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=cellranger_ref_gen_bbb_hg19
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#### SBATCH --mem-per-cpu=30000 # memory limit per core
#SBATCH --mem=50G # memory limit per compute node for the job
#SBATCH --time=3-00:00 # maximum job time in D-HH:MM
#SBATCH --account=scw1329
#SBATCH -o /scratch/c.mpmgb/hawk_output/%x_out_%J.txt
#SBATCH -e /scratch/c.mpmgb/hawk_output/%x_err_%J.txt
#SBATCH --mail-user Bernardo-HarringtonG@cardiff.ac.uk # email on fail
#SBATCH --mail-type END,FAIL

echo "*****************************************************************"
echo "job ID: "$SLURM_JOBID
echo "job name: "$SLURM_JOB_NAME
echo "Run on host: "`hostname`
echo "Number of threads (nproc): "`nproc`
echo "Total memory in GB: "`free -g | grep -oP '\d+' | sed -n 1p`
echo "Used memory in GB: "`free -g | grep -oP '\d+' | sed -n 2p`
echo "Free memory in GB: "`free -g | grep -oP '\d+' | sed -n 3p`
echo "Username: "`whoami`
echo "Started at: "`date`
echo -e "*****************************************************************\n"

#-----------------------------------------------------------------------

## code from Cell Ranger (CR) website
## modified to retrieve GRCh38 genome build from specific ensembl release
## the default CR reference from July 2020 (refdata-gex-GRCh38-2020-A.tar.gz) is based on ensembl v98/Gencode v32
## https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build
##
## NOTE: run first download part on login node
## pre-processing and cellranger mkref on compute node

#-----------------------------------------------------------------------

## CR v7.1.0
CELL_RANGER="/scratch/c.mpmgb/tools/cellranger-7.1.0/bin/cellranger"

## which ensembl and gencode version
ENSEMBL_VER="87"
#GENCODE_VER="43"

## store reference data
REF_DIR="/gluster/dri02/rdscw/shared/webber/reference_genomes"
cd $REF_DIR

# Genome metadata
genome="GRCh37"

# Download source files if they do not exist in reference_sources/ folder
source="reference_sources_ensembl_"${ENSEMBL_VER}
mkdir -p "$source"

fasta_url="http://ftp.ensembl.org/pub/grch37/release-"${ENSEMBL_VER}"/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
fasta_in="${source}/Homo_sapiens.GRCh37.dna.primary_assembly.fa"
gtf_url="http://ftp.ensembl.org/pub/grch37/release-"${ENSEMBL_VER}"/gtf/homo_sapiens/Homo_sapiens.GRCh37."${ENSEMBL_VER}".gtf.gz"
gtf_in="${source}/ensembl.v"${ENSEMBL_VER}".primary_assembly.annotation.gtf"

## compute node dir
COM_DIR="/scratch/scw1329/gmbh/c.mpmgb/blood-brain-barrier-in-ad/03_data/992_genotyping_data/02_cellranger_reference"

#-----------------------------------------------------------------------
## run on login node

# if [ ! -f "$fasta_in" ]; then
#     curl -sS "$fasta_url" | zcat > "$fasta_in"
# fi
# if [ ! -f "$gtf_in" ]; then
#     curl -sS "$gtf_url" | zcat > "$gtf_in"
# fi

#-----------------------------------------------------------------------

## copy the downloaded reference 
# cp -a ${REF_DIR}"/." $COM_DIR

## run on compute node

# Move to the scratch dir
cd $COM_DIR

echo "Run Cell Ranger mkref..."  

$CELL_RANGER mkgtf $gtf_in Homo_sapiens.GRCh37.87.filtered.gtf \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense \
                   --memgb=250 \
                   --nthreads=40


$CELL_RANGER mkref --genome=hg19 \
                   --fasta="$fasta_in" \
                   --genes=Homo_sapiens.GRCh37.87.filtered.gtf \
                   --ref-version=3.0.0 \
                   --memgb=250 \
                   --nthreads=40
