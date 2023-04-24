#!/bin/bash

#SBATCH -p c_highmem_dri1 	###c_highmem_dri1 ## dev, compute, htc, highmem
#SBATCH --job-name=cellranger_ref_gen_bbb
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#### SBATCH --mem-per-cpu=30000 # memory limit per core
#SBATCH --mem=260G # memory limit per compute node for the job
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
ENSEMBL_VER="109"
GENCODE_VER="43"

## store reference data
REF_DIR="/gluster/dri02/rdscw/shared/webber/reference_genomes"
cd $REF_DIR

# Genome metadata
genome="GRCh38"
version="2023_"${ENSEMBL_VER}

# Set up source and build directories
## need to update build dir!
build="GRCh38_2023_ensembl"${ENSEMBL_VER}
mkdir -p "$build"

# Download source files if they do not exist in reference_sources/ folder
source="reference_sources_ensembl"${ENSEMBL_VER}
mkdir -p "$source"

fasta_url="http://ftp.ensembl.org/pub/release-"${ENSEMBL_VER}"/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"${GENCODE_VER}"/gencode.v"${GENCODE_VER}".primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.v"${GENCODE_VER}".primary_assembly.annotation.gtf"

## compute node dir
COM_DIR="/scratch/c.mpmgb/blood-brain-barrier-in-ad/03_data/990_processed_data/001_snrnaseq/03_cellranger_reference"

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
# cp -r $REF_DIR $COM_DIR

## run on compute node

# Move to the scratch dir
cd $COM_DIR

# update build dir
build="reference_genomes/"${build}

# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"


# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_modified="$build/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"


# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


echo "Filter GTF..."

# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    >> "$gtf_filtered"
    
    
echo "Run Cell Ranger mkref..."  

## Create reference package
$CELL_RANGER mkref \
--ref-version="$version" \
--genome="$genome" \
--fasta="$fasta_modified" \
--genes="$gtf_filtered" \
--memgb=250 \
--nthreads=40
