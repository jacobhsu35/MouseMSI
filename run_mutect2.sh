#!/bin/bash
#PBS -l select=1:ncpus=10
#PBS -q ngs48G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N SKLS518_20200717_test
#PBS -o /project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/SKLS518_20200717_test/20200721_Mutect2.out
#PBS -e /project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/SKLS518_20200717_test/20200721_Mutect2.err
#PBS -M jacobhsu@ntu.edu.tw
#PBS -m e
GATK_PATH=/pkg/biology/GATK/GATK_v4.1.8.0
PICARD_PATH=/pkg/biology/Picard/Picard_v2.18.11/picard.jar
BWA_PATH=/pkg/biology/BWA/BWA_v0.7.17
SAMTOOLS_PATH=/pkg/biology/SAMtools/SAMtools_v1.10/bin
BCFTOOLS_PATH=/pkg/biology/BCFtools/BCFtools_v1.10.2/bin
SVABA_PATH=/pkg/biology/SvABA/SvABA_v1.1.0/bin

TUMOR_FASTQ_1_PATH=$1
TUMOR_FASTQ_2_PATH=$2
NORMAL_FASTQ_1_PATH=$3
NORMAL_FASTQ_2_PATH=$4
TUMOR_ID=$5
NORMAL_ID=$6
OUTPUT_PATH=$7
REF_GENOME_PATH=$8
HUMAN_DBSNP_PATH=$9
BAIT=${10}
TARGET=${11}
#GERMLINE_RESOURCE_PATH=${10}
#GERMLINE_RESOURCE_FOR_PILEUP_PATH=${11}

# REF_GENOME_PATH=/project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/hg38/Homo_sapiens_assembly38.fasta
# REF_GENOME_PATH=/home/hsiaoyi0504/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# HUMAN_DBSNP_PATH=/home/hsiaoyi0504/Homo_sapiens_assembly38.dbsnp138.vcf

# Ideally, you don't need to modify following lines

NUM_THREAD=10

mkdir -p $OUTPUT_PATH
cd $OUTPUT_PATH
if [[ $NORMAL_FASTQ_1_PATH == *.fastq.gz ]] || [[ $NORMAL_FASTQ_1_PATH == *.fastq ]]
then
    # for Tumor samples
    $BWA_PATH/bwa mem -Y \
    -R "@RG\tID:GP-${TUMOR_ID}\tSM:${TUMOR_ID}\tPL:illumina" \
    -t $NUM_THREAD \
    -K 100000000 \
    $REF_GENOME_PATH \
    $TUMOR_FASTQ_1_PATH \
    $TUMOR_FASTQ_2_PATH 2> $OUTPUT_PATH/${TUMOR_ID}_bwa.log \
    | $SAMTOOLS_PATH/samtools sort -@$NUM_THREAD -o $OUTPUT_PATH/${TUMOR_ID}.bam

    # for Normal samples
    $BWA_PATH/bwa mem -Y \
    	-R "@RG\tID:GP-${NORMAL_ID}\tSM:${NORMAL_ID}\tPL:illumina" \
    	-t $NUM_THREAD \
        -K 100000000 \
        $REF_GENOME_PATH \
        $NORMAL_FASTQ_1_PATH \
        $NORMAL_FASTQ_2_PATH 2> $OUTPUT_PATH/${NORMAL_ID}_bwa.log \
        | $SAMTOOLS_PATH/samtools sort -@$NUM_THREAD -o $OUTPUT_PATH/${NORMAL_ID}.bam
elif [[ $NORMAL_FASTQ_1_PATH == *.fastq.bz2 ]]
then
    # deal with bz2 files
    # http://seqanswers.com/forums/showthread.php?t=42300

    # for Tumor samples
    $BWA_PATH/bwa mem -Y \
        -R "@RG\tID:GP-${TUMOR_ID}\tSM:${TUMOR_ID}\tPL:illumina" \
        -t $NUM_THREAD \
        -K 100000000 \
        $REF_GENOME_PATH \
        <(bunzip2 -c $TUMOR_FASTQ_1_PATH) \
        <(bunzip2 -c $TUMOR_FASTQ_2_PATH) 2> $OUTPUT_PATH/${TUMOR_ID}_bwa.log \
        | $SAMTOOLS_PATH/samtools sort -@$NUM_THREAD -o $OUTPUT_PATH/${TUMOR_ID}.bam

    # for Normal samples
    $BWA_PATH/bwa mem -Y \
        -R "@RG\tID:GP-${NORMAL_ID}\tSM:${NORMAL_ID}\tPL:illumina" \
        -t $NUM_THREAD \
        -K 100000000 \
        $REF_GENOME_PATH \
        <(bunzip2 -c $NORMAL_FASTQ_1_PATH) \
        <(bunzip2 -c $NORMAL_FASTQ_2_PATH) 2> $OUTPUT_PATH/${NORMAL_ID}_bwa.log \
        | $SAMTOOLS_PATH/samtools sort -@$NUM_THREAD -o $OUTPUT_PATH/${NORMAL_ID}.bam

fi

$GATK_PATH/gatk --java-options "-Xmx40G" ValidateSamFile \
 -I $OUTPUT_PATH/${NORMAL_ID}.bam \
 -O $OUTPUT_PATH/${NORMAL_ID}_validate_sam.log \
 -R $REF_GENOME_PATH \
 -M SUMMARY

$GATK_PATH/gatk --java-options "-Xmx40G" ValidateSamFile \
 -I $OUTPUT_PATH/${TUMOR_ID}.bam \
 -O $OUTPUT_PATH/${TUMOR_ID}_validate_sam.log \
 -R $REF_GENOME_PATH \
 -M SUMMARY 

### mark duplicate (for normal sample)
$GATK_PATH/gatk --java-options "-Xmx40G" MarkDuplicates \
 -I $OUTPUT_PATH/${NORMAL_ID}.bam \
 -O $OUTPUT_PATH/${NORMAL_ID}_marked.bam \
 --METRICS_FILE $OUTPUT_PATH/${NORMAL_ID}_MarkDuplicates_metrics.txt > $OUTPUT_PATH/MarkDuplicates.${TUMOR_ID}.log 2>&1

 ### mark duplicate (for tumor sample)
$GATK_PATH/gatk --java-options "-Xmx40G" MarkDuplicates \
 -I $OUTPUT_PATH/${TUMOR_ID}.bam \
 -O $OUTPUT_PATH/${TUMOR_ID}_marked.bam \
 --METRICS_FILE $OUTPUT_PATH/${TUMOR_ID}_MarkDuplicates_metrics.txt >> $OUTPUT_PATH/MarkDuplicates.${TUMOR_ID}.log 2>&1

### Base Quality Score Recalibration (BQSR) score

### for normal
### BQSR first pass
$GATK_PATH/gatk --java-options "-Xmx40G" BaseRecalibrator \
 -I $OUTPUT_PATH/${NORMAL_ID}_marked.bam \
 -R $REF_GENOME_PATH \
 --known-sites $HUMAN_DBSNP_PATH \
 -O $OUTPUT_PATH/${NORMAL_ID}_recal_pass1.table > $OUTPUT_PATH/BQSR.${TUMOR_ID}.log 2>&1

$GATK_PATH/gatk --java-options "-Xmx40G" ApplyBQSR \
 -I $OUTPUT_PATH/${NORMAL_ID}_marked.bam \
 -R $REF_GENOME_PATH \
 --bqsr-recal-file $OUTPUT_PATH/${NORMAL_ID}_recal_pass1.table \
 -O $OUTPUT_PATH/${NORMAL_ID}_marked.recal.pass1.bam >> $OUTPUT_PATH/BQSR.${TUMOR_ID}.log 2>&1

### BQSR second pass
$GATK_PATH/gatk --java-options "-Xmx40G" BaseRecalibrator \
 -I $OUTPUT_PATH/${NORMAL_ID}_marked.recal.pass1.bam \
 -R $REF_GENOME_PATH \
 --known-sites $HUMAN_DBSNP_PATH \
 -O $OUTPUT_PATH/${NORMAL_ID}_recal_pass2.table >> $OUTPUT_PATH/BQSR.${TUMOR_ID}.log 2>&1

### Analyze covariates
$GATK_PATH/gatk --java-options "-Xmx40G" AnalyzeCovariates \
 --before-report-file $OUTPUT_PATH/${NORMAL_ID}_recal_pass1.table \
 --after-report-file $OUTPUT_PATH/${NORMAL_ID}_recal_pass2.table \
 --plots-report-file $OUTPUT_PATH/${NORMAL_ID}_covariates.pdf > $OUTPUT_PATH/Analyze_covariates.log 2>&1

### for tumor
### BQSR first pass
$GATK_PATH/gatk --java-options "-Xmx40G" BaseRecalibrator \
 -I $OUTPUT_PATH/${TUMOR_ID}_marked.bam \
 -R $REF_GENOME_PATH \
 --known-sites $HUMAN_DBSNP_PATH \
 -O $OUTPUT_PATH/${TUMOR_ID}_recal_pass1.table >> $OUTPUT_PATH/BQSR.${TUMOR_ID}.log 2>&1

$GATK_PATH/gatk --java-options "-Xmx40G" ApplyBQSR \
 -I $OUTPUT_PATH/${TUMOR_ID}_marked.bam \
 -R $REF_GENOME_PATH \
 --bqsr-recal-file $OUTPUT_PATH/${TUMOR_ID}_recal_pass1.table \
 -O $OUTPUT_PATH/${TUMOR_ID}_marked.recal.pass1.bam >> $OUTPUT_PATH/BQSR.${TUMOR_ID}.log 2>&1

 ### BQSR second pass
$GATK_PATH/gatk --java-options "-Xmx40G" BaseRecalibrator \
 -I $OUTPUT_PATH/${TUMOR_ID}_marked.recal.pass1.bam \
 -R $REF_GENOME_PATH \
 --known-sites $HUMAN_DBSNP_PATH \
 -O $OUTPUT_PATH/${TUMOR_ID}_recal_pass2.table >> $OUTPUT_PATH/BQSR.${TUMOR_ID}.log 2>&1

 ### Analyze covariates
$GATK_PATH/gatk --java-options "-Xmx40G" AnalyzeCovariates \
 --before-report-file $OUTPUT_PATH/${TUMOR_ID}_recal_pass1.table \
 --after-report-file $OUTPUT_PATH/${TUMOR_ID}_recal_pass2.table \
 --plots-report-file $OUTPUT_PATH/${TUMOR_ID}_covariates.pdf >> $OUTPUT_PATH/Analyze_covariates.log 2>&1

### Mutect
### get sample names that will be used for later in several command line calls 

$GATK_PATH/gatk --java-options "-Xmx40G" GetSampleName -I $OUTPUT_PATH/${TUMOR_ID}_marked.recal.pass1.bam -O $OUTPUT_PATH/${TUMOR_ID}_sample_name.txt
$GATK_PATH/gatk --java-options "-Xmx40G" GetSampleName -I $OUTPUT_PATH/${NORMAL_ID}_marked.recal.pass1.bam -O $OUTPUT_PATH/${NORMAL_ID}_sample_name.txt

### Calculate Metrics of inverval

TUMOR_SAMPLE_NAME=$(cat $OUTPUT_PATH/${TUMOR_ID}_sample_name.txt)
NORMAL_SAMPLE_NAME=$(cat $OUTPUT_PATH/${NORMAL_ID}_sample_name.txt)

java -jar $PICARD_PATH CollectHsMetrics I=$OUTPUT_PATH/${TUMOR_ID}_marked.recal.pass1.bam o=$OUTPUT_PATH/${TUMOR_SAMPLE_NAME}.metrics.txt R=${REF_GENOME_PATH} BAIT_INTERVALS=${BAIT} TARGET_INTERVALS=${TARGET} > $TUMOR_SAMPLE_NAME"_"collecthsmetrics.log 2>&1
java -jar $PICARD_PATH CollectHsMetrics I=$OUTPUT_PATH/${NORMAL_ID}_marked.recal.pass1.bam o=$OUTPUT_PATH/${NORMAL_SAMPLE_NAME}.metrics.txt R=${REF_GENOME_PATH} BAIT_INTERVALS=${BAIT} TARGET_INTERVALS=${TARGET} > $NORMAL_SAMPLE_NAME"_"collecthsmetrics.log 2>&1

### Call somatic short variants and generate a BAM with Mutect2

$GATK_PATH/gatk --java-options "-Xmx40G" Mutect2 \
 -R $REF_GENOME_PATH \
 -I $OUTPUT_PATH/${TUMOR_ID}_marked.recal.pass1.bam \
 -I $OUTPUT_PATH/${NORMAL_ID}_marked.recal.pass1.bam \
 -tumor $TUMOR_SAMPLE_NAME \
 -normal $NORMAL_SAMPLE_NAME \
 -O $OUTPUT_PATH/Mutect2.vcf.gz \
 -minimum-allele-fraction 0.1 \
 -enable-all-annotations true \
 --annotations-to-exclude UniqueAltReadCount \
 --f1r2-tar-gz $OUTPUT_PATH/f1r2.tar.gz \
 --native-pair-hmm-threads $NUM_THREAD \
 --create-output-bam-md5 true \
 -L ${TARGET} \
 -bamout $OUTPUT_PATH/Mutect2.bam > $OUTPUT_PATH/mutect2.log 2>&1
# -germline-resource $GERMLINE_RESOURCE_PATH \  https://gatk.broadinstitute.org/hc/en-us/community/posts/360056904972-JAVA-errors-from-gnomad-resource-in-Mutect2

#$GATK_PATH/gatk --java-options "-Xmx40G" LearnReadOrientationModel -I $OUTPUT_PATH/f1r2.tar.gz -O $OUTPUT_PATH/read-orientation-model.tar.gz

#$GATK_PATH/gatk --java-options "-Xmx40G" GetPileupSummaries \
#    -I $OUTPUT_PATH/${TUMOR_ID}_marked.recal.pass1.bam \
#    -V $GERMLINE_RESOURCE_FOR_PILEUP_PATH \
#    -L $GERMLINE_RESOURCE_FOR_PILEUP_PATH \
#    -O $OUTPUT_PATH/getpileupsummaries.table > $OUTPUT_PATH/get_pileup_summaries.log 2>&1

#$GATK_PATH/gatk --java-options "-Xmx40G" CalculateContamination \
#    -I $OUTPUT_PATH/getpileupsummaries.table \
#    -tumor-segmentation $OUTPUT_PATH/${TUMOR_ID}_segments.table \
#    -O $OUTPUT_PATH/contamination.table > $OUTPUT_PATH/calculate_contamination.log 2>&1

$GATK_PATH/gatk --java-options "-Xmx40G" FilterMutectCalls \
    -R $REF_GENOME_PATH \
    -V $OUTPUT_PATH/Mutect2.vcf.gz \
    -O $OUTPUT_PATH/Mutect2.filtered.vcf.gz \
    --min-allele-fraction 0.01 \
    --unique-alt-read-count 10 > $OUTPUT_PATH/filter_mutect_calls.log 2>&1
#    --tumor-segmentation $OUTPUT_PATH/${TUMOR_ID}_segments.table \
#    --contamination-table $OUTPUT_PATH/contamination.table \
#    --ob-priors $OUTPUT_PATH/read-orientation-model.tar.gz \

gzip -d $OUTPUT_PATH/Mutect2.filtered.vcf.gz
VT_PATH=/project/GP1/u3710062/AI_SHARE/software/vt-0.57721
$VT_PATH/vt decompose -s -o Mutect2.decom.vcf Mutect2.filtered.vcf
$VT_PATH/vt normalize -o Mutect2.norm.vcf -r $REF_GENOME_PATH Mutect2.decom.vcf

$GATK_PATH/gatk --java-options "-Xmx40G" SelectVariants \
        -R $REF_GENOME_PATH \
        -V $OUTPUT_PATH/Mutect2.norm.vcf \
        --select-type-to-include SNP \
        -O $OUTPUT_PATH/Mutect2.norm.snp.vcf.gz

$GATK_PATH/gatk --java-options "-Xmx40G" SelectVariants \
        -R $REF_GENOME_PATH \
        -V $OUTPUT_PATH/Mutect2.norm.vcf \
        --select-type-to-include INDEL \
        -O $OUTPUT_PATH/Mutect2.norm.indel.vcf.gz

$GATK_PATH/gatk --java-options "-Xmx40G" SelectVariants \
        -R $REF_GENOME_PATH \
        -V $OUTPUT_PATH/Mutect2.norm.vcf \
        --exclude-filtered true \
        -O $OUTPUT_PATH/Mutect2.norm.PASS.vcf.gz

#cd $OUTPUT_PATH

#$SVABA_PATH/svaba run -t $OUTPUT_PATH/${TUMOR_ID}_marked.recal.pass1.bam -n $OUTPUT_PATH/${NORMAL_ID}_marked.recal.pass1.bam -G $REF_GENOME_PATH -a somatic_run -p $NUM_THREAD -D /project/GP1/u3710062/AI_SHARE/reference/GRCm38_dbSNP/00-All.vcf.gz 

#SVABA_INDEL_VCF_PATH=$OUTPUT_PATH/somatic_run.svaba.somatic.indel.vcf

#$BCFTOOLS_PATH/bcftools view $SVABA_INDEL_VCF_PATH -Oz -o $OUTPUT_PATH/svaba.indel.vcf.gz
#$BCFTOOLS_PATH/bcftools index $OUTPUT_PATH/svaba.indel.vcf.gz

#$BCFTOOLS_PATH/bcftools isec -n=2 $OUTPUT_PATH/Mutect2.indel.vcf.gz $OUTPUT_PATH/svaba.indel.vcf.gz -w 1 -Oz -o $OUTPUT_PATH/consensus.indel.vcf.gz

