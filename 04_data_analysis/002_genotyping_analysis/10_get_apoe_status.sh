#!/bin/bash

# Set the path to your VCF files directory
vcf_dir="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/992_genotyping_data/PLINK_021123_0956/split_samples"

# Set the output file path
output_file="/scratch/scw1329/gmbh/blood-brain-barrier-in-ad/03_data/992_genotyping_data/apoe_status.csv"

# Write the CSV header
echo "sample_name,apoe_status,rs7412,rs429358" > "$output_file"

# Loop through each VCF file in the directory
for vcf_file in "$vcf_dir"/*.vcf; do
    # Extract sample name from the file name
    sample_name=$(basename "$vcf_file" .vcf)

    # Use awk to extract APOE rs7412 and rs429358 genotypes for the sample
    # note that cases where the ALT is "." are excluded
    apoE_rs7412=$(awk '$1 == "chr19" && $3 ~ /rs7412/ && $5 != "." {print $10}' "$vcf_file")
    apoE_rs429358=$(awk '$1 == "chr19" && $3 ~ /rs429358/ && $5 != "." {print $10}' "$vcf_file")


    # Count occurrences of each genotype
    count_rs7412_0_0=$(echo "$apoE_rs7412" | grep -o '0/0' | wc -l)
    count_rs7412_0_1=$(echo "$apoE_rs7412" | grep -o '0/1' | wc -l)
    count_rs7412_1_1=$(echo "$apoE_rs7412" | grep -o '1/1' | wc -l)
    count_rs7412_1_0=$(echo "$apoE_rs7412" | grep -o '1/0' | wc -l)

    count_rs429358_0_0=$(echo "$apoE_rs429358" | grep -o '0/0' | wc -l)
    count_rs429358_0_1=$(echo "$apoE_rs429358" | grep -o '0/1' | wc -l)
    count_rs429358_1_1=$(echo "$apoE_rs429358" | grep -o '1/1' | wc -l)
    count_rs429358_1_0=$(echo "$apoE_rs429358" | grep -o '1/0' | wc -l)

    # Determine the majority genotype
    majority_genotype_rs7412=$(echo "$count_rs7412_0_0 $count_rs7412_0_1 $count_rs7412_1_1 $count_rs7412_1_0" | awk '{max=$1;genotype="0/0";for(i=2;i<=NF;i++)if($i>max){max=$i;genotype="0/0"};print genotype}')
    majority_genotype_rs429358=$(echo "$count_rs429358_0_0 $count_rs429358_0_1 $count_rs429358_1_1 $count_rs429358_1_0" | awk '{max=$1;genotype="0/0";for(i=2;i<=NF;i++)if($i>max){max=$i;genotype="0/0"};print genotype}')

    #Translate genotypes to APOE status
    case "$majority_genotype_rs7412$majority_genotype_rs429358" in
        "1/10/0") apoE_status="E2/E2" ;;
        "0/00/0") apoE_status="E3/E3" ;;
        "1/00/0") apoE_status="E2/E3" ;;
        "0/10/0") apoE_status="E3/E2" ;;
        "0/11/0") apoE_status="E4/E2" ;;
        "1/00/1") apoE_status="E2/E4" ;;
        "0/01/0") apoE_status="E4/E3" ;;
        "0/01/1") apoE_status="E4/E4" ;;
        "0/00/1") apoE_status="E3/E4" ;;
        *) apoE_status="Unknown" ;;
    esac

    # Echo the majority allele for each SNP
    echo "Sample: $sample_name, APOE rs7412: $majority_genotype_rs7412, APOE rs429358: $majority_genotype_rs429358"
    # Output the APOE status for the sample
    echo "$sample_name,$apoE_status,$majority_genotype_rs7412,$majority_genotype_rs429358" >> "$output_file"
done

echo "Results saved to $output_file"
