#!/bin/bash


## collect MultiQC reports
OUTPUT_DIR="/Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/MultiQC_reports/"


#-----------------------------------------------------------------------
# FastQC runs

## Endo 10X Vascular and Parenchymal fraction set 1 - 24 samples
## INPUT_DIR="/Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/FastQC/FastQC_main_set1/"


#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_I1" \
#--filename "fastqc_endo_10X_set1_I1" \
#--title "FastQC endo 10X set1 I1" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""*_I1_*_fastqc.zip

#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_R1" \
#--filename "fastqc_endo_10X_set1_R1" \
#--title "FastQC endo 10X set1 R1" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""*_R1_*_fastqc.zip

#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_R2" \
#--filename "fastqc_endo_10X_set1_R2" \
#--title "FastQC endo 10X set1 R2" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""*_R2_*_fastqc.zip


## combine all 4 runs from set 2 in one report
#INPUT_DIR="/Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/FastQC/Set_2/"

#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set2_I1" \
#--filename "fastqc_endo_10X_set2_I1" \
#--title "FastQC endo 10X set2 I1" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""/*/*_I1_*_fastqc.zip


#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set2_I2" \
#--filename "fastqc_endo_10X_set2_I2" \
#--title "FastQC endo 10X set2 I2" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""/*/*_I2_*_fastqc.zip

#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set2_R1" \
#--filename "fastqc_endo_10X_set2_R1" \
#--title "FastQC endo 10X set2 R1" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""/*/*_R1_*_fastqc.zip


#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set2_R2" \
#--filename "fastqc_endo_10X_set2_R2" \
#--title "FastQC endo 10X set2 R2" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""/*/*_R2_*_fastqc.zip


#INPUT_DIR="/Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/FastQC/FastQC_main_set1_miseq/"

#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_R1_miseq" \
#--filename "fastqc_endo_10X_set1_R1_miseq" \
#--title "FastQC endo 10X set1 R1" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""*_R1_*_fastqc.zip

#multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_R2_miseq" \
#--filename "fastqc_endo_10X_set1_R2)miseq" \
#--title "FastQC endo 10X set1 R2" \
#--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
#--module fastqc \
#--force \
#--interactive \
#$INPUT_DIR""*_R2_*_fastqc.zip




INPUT_DIR="/Users/frankwessely/Desktop/Documents/Projects/BBB/Endo_10X_Main/Analysis/FastQC/FastQC_main_set1_v2/"

multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_v2_I1" \
--filename "fastqc_endo_10X_set1_I1" \
--title "FastQC endo 10X set1 v2 I1" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR""*_I1_*_fastqc.zip


multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_v2_I2" \
--filename "fastqc_endo_10X_set1_I2" \
--title "FastQC endo 10X set1 v2 I2" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR""*_I2_*_fastqc.zip

multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_v2_R1" \
--filename "fastqc_endo_10X_set1_R1" \
--title "FastQC endo 10X set1 v2 R1" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR""*_R1_*_fastqc.zip


multiqc -o $OUTPUT_DIR"FastQC_endo_10X_set1_v2_R2" \
--filename "fastqc_endo_10X_set1_R2" \
--title "FastQC endo 10X set1 v2 R2" \
--cl_config "fastqc_config: { fastqc_theoretical_gc: hg38_txome }" \
--module fastqc \
--force \
--interactive \
$INPUT_DIR""*_R2_*_fastqc.zip

