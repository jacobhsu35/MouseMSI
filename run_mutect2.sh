#!/bin/bash

GATK_PATH=/pkg/biology/GATK/GATK_v4.1.6.0
PICARD_PATH=/pkg/biology/Picard/Picard_v2.18.11/picard.jar
BWA_PATH=/pkg/biology/BWA/BWA_v0.7.17
SAMTOOLS_PATH=/pkg/biology/SAMtools/SAMtools_v1.10/bin

REF_GENOME_PATH=/project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/hg38/Homo_sapiens_assembly38.fasta
# REF_GENOME_PATH=/home/hsiaoyi0504/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
HUMAN_DBSNP_PATH=/home/hsiaoyi0504/Homo_sapiens_assembly38.dbsnp138.vcf

TUMOR_FASTQ_1_PATH=$1
TUMOR_FASTQ_2_PATH=$2
NORMAL_FASTQ_1_PATH=$3
NORMAL_FASTQ_2_PATH=$4
TUMOR_READGROUP=$5
NORMAL_READGROUP=$6
OUTPUT_PATH=$7


# Ideally, you don't need to modify following lines

NUM_THREAD=10

mkdir -p $OUTPUT_PATH

if [[ $NORMAL_FASTQ_1_PATH == *.fastq.gz ]] || [[ $NORMAL_FASTQ_1_PATH == *.fastq ]]
then
    # for normal samples
    $BWA_PATH/bwa mem -M \
        -R $NORMAL_READGROUP \
        -t $NUM_THREAD \
        -K 10000000 \
        $REF_GENOME_PATH \
        $NORMAL_FASTQ_1_PATH \
        $NORMAL_FASTQ_2_PATH 2> $OUTPUT_PATH/normal_bwa.log \
        | $SAMTOOLS_PATH/samtools sort -@10 -o $OUTPUT_PATH/normal.bam

    # for tumor samples
    $BWA_PATH/bwa mem -M \
        -R $TUMOR_READGROUP \
        -t $NUM_THREAD \
        -K 10000000 \
        $REF_GENOME_PATH \
        $TUMOR_FASTQ_1_PATH \
        $TUMOR_FASTQ_2_PATH 2> $OUTPUT_PATH/tumor_bwa.log \
        | $SAMTOOLS_PATH/samtools sort -@10 -o $OUTPUT_PATH/tumor.bam
elif [[ $NORMAL_FASTQ_1_PATH == *.fastq.bz2 ]]
then
    # deal with bz2 files
    # http://seqanswers.com/forums/showthread.php?t=42300
    # for normal samples
    $BWA_PATH/bwa mem -M \
        -R $NORMAL_READGROUP \
        -t $NUM_THREAD \
        -K 10000000 \
        $REF_GENOME_PATH \
        <(bunzip2 -c $NORMAL_FASTQ_1_PATH) \
        <(bunzip2 -c $NORMAL_FASTQ_2_PATH) 2> $OUTPUT_PATH/normal_bwa.log \
        | $SAMTOOLS_PATH/samtools sort -@10 -o $OUTPUT_PATH/normal.bam

    # for tumor samples
    $BWA_PATH/bwa mem -M \
        -R $TUMOR_READGROUP \
        -t $NUM_THREAD \
        -K 10000000 \
        $REF_GENOME_PATH \
        <(bunzip2 -c $TUMOR_FASTQ_1_PATH) \
        <(bunzip2 -c $TUMOR_FASTQ_2_PATH) 2> $OUTPUT_PATH/tumor_bwa.log \
        | $SAMTOOLS_PATH/samtools sort -@10 -o $OUTPUT_PATH/tumor.bam
fi

java -jar $PICARD_PATH ValidateSamFile \
INPUT=$OUTPUT_PATH/normal.bam \
OUTPUT=$OUTPUT_PATH/normal_validate_sam.log \
REFERENCE_SEQUENCE=$REF_GENOME_PATH \
MODE=SUMMARY

java -jar $PICARD_PATH ValidateSamFile \
INPUT=$OUTPUT_PATH/tumor.bam \
OUTPUT=$OUTPUT_PATH/tumor_validate_sam.log \
REFERENCE_SEQUENCE=$REF_GENOME_PATH \
MODE=SUMMARY 

# mark duplicate (for normal sample)
java -jar $PICARD_PATH MarkDuplicates \
 INPUT=$OUTPUT_PATH/normal.bam \
 OUTPUT=$OUTPUT_PATH/normal_marked.bam \
 METRICS_FILE=$OUTPUT_PATH/normal_metrics.txt > $OUTPUT_PATH/picard_normal.log 2>&1

# mark duplicate (for tumor sample)
java -jar $PICARD_PATH MarkDuplicates \
 INPUT=$OUTPUT_PATH/tumor.bam \
 OUTPUT=$OUTPUT_PATH/tumor_marked.bam \
 METRICS_FILE=$OUTPUT_PATH/tumor_metrics.txt > $OUTPUT_PATH/picard_tumor.log 2>&1

# Base Quality Score Recalibration (BQSR) score

## for normal
### BQSR first pass
$GATK_PATH/gatk BaseRecalibrator \
 -I=$OUTPUT_PATH/normal_marked.bam \
 -R=$REF_GENOME_PATH \
 --known-sites $HUMAN_DBSNP_PATH \
 -O=$OUTPUT_PATH/normal_recal_pass1.table > $OUTPUT_PATH/normal_BQSR_first_pass.log 2>&1

$GATK_PATH/gatk ApplyBQSR \
 -I=$OUTPUT_PATH/normal_marked.bam \
 -R=$REF_GENOME_PATH \
 --bqsr-recal-file $OUTPUT_PATH/normal_recal_pass1.table \
 -O=$OUTPUT_PATH/normal_marked.recal.pass1.bam > $OUTPUT_PATH/normal_apply_BQSR.log 2>&1

### BQSR second pass
$GATK_PATH/gatk BaseRecalibrator \
 -I=$OUTPUT_PATH/normal_marked.recal.pass1.bam \
 -R=$REF_GENOME_PATH \
 --known-sites $HUMAN_DBSNP_PATH \
 -O=$OUTPUT_PATH/normal_recal_pass2.table > $OUTPUT_PATH/normal_BQSR_second_pass.log 2>&1

### Analyze covariates
$GATK_PATH/gatk AnalyzeCovariates \
 --before-report-file $OUTPUT_PATH/normal_recal_pass1.table \
 --after-report-file $OUTPUT_PATH/normal_recal_pass2.table \
 --plots-report-file $OUTPUT_PATH/normal_covariates.pdf 

## for tumor
### BQSR first pass
$GATK_PATH/gatk BaseRecalibrator \
 -I=$OUTPUT_PATH/tumor_marked.bam \
 -R=$REF_GENOME_PATH \
 --known-sites $HUMAN_DBSNP_PATH \
 -O=$OUTPUT_PATH/tumor_recal_pass1.table

$GATK_PATH/gatk ApplyBQSR \
 -I=$OUTPUT_PATH/tumor_marked.bam \
 -R=$REF_GENOME_PATH \
 --bqsr-recal-file $OUTPUT_PATH/tumor_recal_pass1.table \
 -O=$OUTPUT_PATH/tumor_marked.recal.pass1.bam

 ### BQSR second pass
$GATK_PATH/gatk BaseRecalibrator \
 -I=$OUTPUT_PATH/tumor_marked.recal.pass1.bam \
 -R=$REF_GENOME_PATH \
 --known-sites $HUMAN_DBSNP_PATH \
 -O=$OUTPUT_PATH/tumor_recal_pass2.table

### Analyze covariates
$GATK_PATH/gatk AnalyzeCovariates \
 --before-report-file $OUTPUT_PATH/tumor_recal_pass1.table \
 --after-report-file $OUTPUT_PATH/tumor_recal_pass2.table \
 --plots-report-file $OUTPUT_PATH/tumor_covariates.pdf

# Mutect
## get sample names that will be used for later in several command line calls 

$GATK_PATH/gatk GetSampleName  -I $OUTPUT_PATH/tumor_marked.recal.pass1.bam -O $OUTPUT_PATH/tumor_sample_name.txt
$GATK_PATH/gatk GetSampleName  -I $OUTPUT_PATH/normal_marked.recal.pass1.bam -O $OUTPUT_PATH/normal_sample_name.txt

TUMOR_SAMPLE_NAME=$(cat $OUTPUT_PATH/tumor_sample_name.txt)
NORMAL_SAMPLE_NAME=$(cat $OUTPUT_PATH/normal_sample_name.txt)

## Call somatic short variants and generate a BAM with Mutect2

$GATK_PATH/gatk Mutect2 \
 -R $REF_GENOME_PATH \
 -I $OUTPUT_PATH/tumor_marked.recal.pass1.bam \
 -I $OUTPUT_PATH/normal_marked.recal.pass1.bam \
 -tumor $TUMOR_SAMPLE_NAME \
 -normal $NORMAL_SAMPLE_NAME \
 -O $OUTPUT_PATH/Mutect2.vcf.gz \
 -bamout $OUTPUT_PATH/Mutect2.bam

