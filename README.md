# Germline Variant Calling Pipeline

## Overview
This pipeline processes paired-end FASTQ files for **germline variant calling** using **BWA-MEM, GATK HaplotypeCaller, and bcftools**. It includes quality control, read alignment, duplicate marking, base recalibration, variant calling, joint genotyping, and variant quality score recalibration (VQSR) to produce high-confidence variant calls.

## Features
- **End-to-end automation** for germline variant calling
- **Quality control and filtering** for improved variant accuracy
- **Joint genotyping** for multiple samples
- **Variant Quality Score Recalibration (VQSR)** for high-confidence calls
- **Logging and error handling** to track each step

## Pipeline Steps

1. **Quality Filtering**: Trimming and filtering of low-quality reads using `fastp`.
2. **Read Alignment**: Mapping reads to the reference genome (`hg38`) using `bwa mem`.
3. **Add Read Groups**: Necessary for GATK compatibility using `AddOrReplaceReadGroups`.
4. **Duplicate Marking**: Identifying duplicate reads with `MarkDuplicatesSpark`.
5. **Base Quality Score Recalibration (BQSR)**: Two-step recalibration using known variant sites.
6. **Variant Calling**: Using `HaplotypeCaller` to generate `.g.vcf.gz` files.
7. **Combine GVCFs**: Merging individual GVCFs into a cohort GVCF.
8. **Joint Genotyping**: Generating a multi-sample VCF.
9. **VQSR (SNP & INDEL)**: Calibrating variant quality for SNPs and INDELs.
10. **VCF Splitting**: Creating individual VCFs from the final cohort VCF.

## Installation
Ensure the following tools are installed:
- **BWA**
- **Samtools**
- **Fastp**
- **GATK**
- **bcftools**

Install example:
```sh
sudo apt update && sudo apt install -y bwa samtools fastp bcftools
# Download GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip
unzip gatk-4.2.6.1.zip && mv gatk-4.2.6.1 /opt/gatk
export PATH=/opt/gatk:$PATH
```

## Usage
```sh
bash germline_variant_calling.sh /path/to/fastq_files
```

## Output Files
- BAM files with alignment statistics
- GVCF files (intermediate)
- Raw and recalibrated cohort VCFs
- Final filtered per-sample VCFs
- Quality control HTML reports
- Log files

## Script
The following is the **full pipeline script**:

```bash
#!/bin/bash
# Germline Variant Calling Pipeline

# Allow the script to continue even if some samples fail
set +e
# Ensure pipeline failures are detected
set -o pipefail

# Set up input and output directories
inputDir=$(pwd)
dir_qc="$inputDir/fastqc"
dir_fastp="$inputDir/filtered_qc_report"
dir_map="$inputDir/Mapsam"
dir_vcf="$inputDir/Germline_VCF"
dir_gvcf="$inputDir/GVCFs"
mkdir -p "$dir_map" "$dir_vcf" "$dir_gvcf"

# Reference genome and known sites
file_ref="/path/to/reference/hg38.fa"
known_site="/path/to/vcf/dbsnp.vcf.gz"
platform="ILLUMINA"

# Logging setup
timestamp=$(date +'%Y%m%d%H%M%S')
log_file="$inputDir/log_${timestamp}.txt"
terminal_log="$inputDir/terminal_log_${timestamp}.txt"
exec > >(tee -a "$terminal_log") 2>&1

echo "Script started." | tee -a "$log_file"

# Detect unique sample prefixes from FASTQ filenames
samples=()
for f in "$inputDir"/*.fastq.gz; do
    if [ -f "$f" ]; then
        sampleName=$(basename "$f" | awk -F"_" '{print $1;}')
        samples+=("$sampleName")
    fi
done

# Remove duplicates to get unique samples
unique_samples=($(printf "%s\n" "${samples[@]}" | sort -u))
numSamples=${#unique_samples[@]}
echo "Total samples: $numSamples" | tee -a "$log_file"

# Function to process each sample through the pipeline
process_sample() {
    local sample="$1"
    echo "Processing sample: $sample" | tee -a "$log_file"

    local filep1="$inputDir/${sample}_1.fastq.gz"
    local filep2="$inputDir/${sample}_2.fastq.gz"

    # Step 1: Perform quality filtering and trimming using Fastp
    fastp -i "$filep1" -I "$filep2" --disable_length_filtering --qualified_quality_phred 30 \
        -o "${dir_fastp}/${sample}_1_filtered.fastq.gz" -O "${dir_fastp}/${sample}_2_filtered.fastq.gz" \
        --html "${dir_fastp}/${sample}.html" --thread 2
    echo "Fastp completed for $sample." | tee -a "$log_file"

    # Step 2: Align reads to reference genome using BWA-MEM
    bwa mem -v 2 -M -t 32 -R "@RG\tID:${sample}\tPL:$platform\tLB:LIB_${sample}\tSM:${sample}" \
        $file_ref ${dir_fastp}/${sample}_1_filtered.fastq.gz ${dir_fastp}/${sample}_2_filtered.fastq.gz | \
        samtools sort -@ 12 -O BAM -o ${dir_map}/${sample}.bam
    echo "BWA MEM completed for $sample." | tee -a "$log_file"

    # Step 3: Generate alignment summary statistics
    samtools flagstat ${dir_map}/${sample}.bam > ${dir_map}/${sample}.Stat.txt

    # Step 4: Add read group information to BAM file
    gatk AddOrReplaceReadGroups -I ${dir_map}/${sample}.bam -O ${dir_map}/${sample}_RG.bam \
        --RGLB Lib-${sample} --RGPL $platform --RGPU ${platform}_${sample} --RGSM $sample \
        --VALIDATION_STRINGENCY LENIENT
    echo "AddOrReplaceReadGroups completed for $sample." | tee -a "$log_file"

    # Step 5: Mark duplicate reads using GATK
    gatk MarkDuplicatesSpark -I ${dir_map}/${sample}_RG.bam -O ${dir_map}/${sample}_markdup.bam \
        --spark-master local[12]
    echo "MarkDuplicatesSpark completed for $sample." | tee -a "$log_file"

    # Step 6: Create base quality score recalibration (BQSR) table
    gatk BaseRecalibrator --java-options "-Xmx4g -XX:ParallelGCThreads=2" \
        -I ${dir_map}/${sample}_markdup.bam --known-sites $known_site \
        -O ${dir_map}/${sample}_recall.table -R $file_ref
    echo "BaseRecalibrator completed for $sample." | tee -a "$log_file"

    # Step 7: Apply BQSR to recalibrate base qualities
    gatk ApplyBQSR --java-options "-Xmx4g -XX:ParallelGCThreads=2" \
        --bqsr-recal-file ${dir_map}/${sample}_recall.table \
        -I ${dir_map}/${sample}_markdup.bam -O ${dir_map}/${sample}_recall.bam -R $file_ref
    echo "ApplyBQSR completed for $sample." | tee -a "$log_file"

    # Step 8: Call variants with HaplotypeCaller in GVCF mode
    gatk HaplotypeCaller --java-options "-Xms20G -Xmx20G -XX:ParallelGCThreads=2" \
        -R $file_ref -I ${dir_map}/${sample}_recall.bam \
        -O ${dir_gvcf}/${sample}.g.vcf.gz -ERC GVCF
    echo "HaplotypeCaller completed for $sample." | tee -a "$log_file"
}

# Step 9: Loop over all samples and call the processing function
for sample in ${unique_samples[@]}; do
    process_sample "$sample"
done

# Step 10: Combine individual GVCFs into one
combined_gvcf="${dir_vcf}/combined.g.vcf.gz"
gatk CombineGVCFs --java-options "-Xms10G -Xmx10G -XX:ParallelGCThreads=4" \
    -R $file_ref $(for sample in ${unique_samples[@]}; do echo "-V ${dir_gvcf}/${sample}.g.vcf.gz"; done) \
    -O $combined_gvcf

# Step 11: Perform joint genotyping on combined GVCFs
cohort_vcf="${dir_vcf}/cohort_raw.vcf.gz"
gatk GenotypeGVCFs --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=4" \
    -R $file_ref -V $combined_gvcf -O $cohort_vcf

# Step 12: SNP recalibration using known training set
gatk VariantRecalibrator --java-options "-Xmx8G -XX:ParallelGCThreads=4" \
    -R $file_ref -V $cohort_vcf --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $training_set_snps \
    -mode SNP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
    -recal-file ${dir_vqsr}/recalibrate_SNP.recal -tranches-file ${dir_vqsr}/recalibrate_SNP.tranches

# Step 13: Apply SNP VQSR to raw cohort VCF
gatk ApplyVQSR --java-options "-Xmx8G -XX:ParallelGCThreads=4" \
    -R $file_ref -V $cohort_vcf --truth-sensitivity-filter-level 99.0 \
    --recal-file ${dir_vqsr}/recalibrate_SNP.recal --tranches-file ${dir_vqsr}/recalibrate_SNP.tranches \
    -mode SNP -O ${dir_vqsr}/cohort_snps.vcf.gz

# Step 14: INDEL recalibration using known indel training set
gatk VariantRecalibrator --java-options "-Xmx8G -XX:ParallelGCThreads=4" \
    -R $file_ref -V ${dir_vqsr}/cohort_snps.vcf.gz \
    --resource:mills,known=true,training=true,truth=true,prior=12.0 $training_set_indels \
    -mode INDEL -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
    -recal-file ${dir_vqsr}/recalibrate_INDEL.recal -tranches-file ${dir_vqsr}/recalibrate_INDEL.tranches

# Step 15: Apply INDEL VQSR to SNP-recalibrated VCF
gatk ApplyVQSR --java-options "-Xmx8G -XX:ParallelGCThreads=4" \
    -R $file_ref -V ${dir_vqsr}/cohort_snps.vcf.gz --truth-sensitivity-filter-level 99.0 \
    --recal-file ${dir_vqsr}/recalibrate_INDEL.recal --tranches-file ${dir_vqsr}/recalibrate_INDEL.tranches \
    -mode INDEL -O ${dir_vcf}/cohort_vqsr.vcf.gz

# Step 16: Split the final cohort VQSR VCF into individual VCFs
for sample in ${unique_samples[@]}; do
    bcftools view -s $sample -Oz -o ${dir_annotated}/${sample}_filtered.vcf.gz ${dir_vcf}/cohort_vqsr.vcf.gz
    echo "Split VCF completed for $sample." | tee -a "$log_file"
done

# Pipeline completion message
echo "Pipeline completed." | tee -a "$log_file"
```

## License
This pipeline is released under the MIT License.

For any issues or improvements, feel free to contribute!
