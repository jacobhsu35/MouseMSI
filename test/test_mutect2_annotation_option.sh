GATK_PATH=/volume/cancer-mutect2/tools/gatk-4.1.5.0
REF_GENOME=hg19_main_chr.fasta
IN_NORMAL_FILE=TM2N.normal.bam
IN_TUMOR_FILE=TM2-2.tumor.bam
GERMLINE_RESOURCE_PATH=af-only-gnomad.raw.sites.b37.germline.vcf.gz
NUM_PAIR_HMM_THREADS=32

# $GATK_PATH/gatk GetSampleName -I $IN_NORMAL_FILE -O tumor_sample_name.txt
# $GATK_PATH/gatk GetSampleName -I $IN_TUMOR_FILE -O normal_sample_name.txt

TUMOR_SAMPLE_NAME=$(cat tumor_sample_name.txt)
NORMAL_SAMPLE_NAME=$(cat normal_sample_name.txt)

$GATK_PATH/gatk Mutect2 \
    -R $REF_GENOME \
    -I $IN_NORMAL_FILE \
    -I $IN_TUMOR_FILE \
    -tumor $TUMOR_SAMPLE_NAME \
    -normal $NORMAL_SAMPLE_NAME \
    -germline-resource $GERMLINE_RESOURCE_PATH \
    -O Mutect2_add_annotation.vcf.gz \
    --f1r2-tar-gz f1r2.tar.gz \
    -bamout Mutect2_result_add_annotation.bam \
    -A AS_BaseQualityRankSumTest \
    -A AS_FisherStrand \
    -A AS_InbreedingCoeff \
    -A AS_MappingQualityRankSumTest \
    -A AS_QualByDepth \
    -A AS_RMSMappingQuality \
    -A AS_ReadPosRankSumTest \
    -A AS_StrandOddsRatio \
    -A AlleleFraction \
    -A BaseQuality \
    -A BaseQualityRankSumTest \
    -A ChromosomeCounts \
    -A ClippingRankSumTest \
    -A CountNs \
    -A Coverage \
    -A DepthPerAlleleBySample \
    -A DepthPerSampleHC \
    -A ExcessHet \
    -A FisherStrand \
    -A FragmentLength \
    -A GenotypeSummaries \
    -A InbreedingCoeff \
    -A LikelihoodRankSumTest \
    -A MappingQuality \
    -A MappingQualityRankSumTest \
    -A MappingQualityZero \
    -A OrientationBiasReadCounts \
    -A OriginalAlignment \
    -A PossibleDeNovo \
    -A QualByDepth \
    -A RMSMappingQuality \
    -A ReadPosRankSumTest \
    -A ReadPosition \
    -A ReferenceBases \
    -A SampleList \
    -A StrandBiasBySample \
    -A StrandOddsRatio \
    -A TandemRepeat \
    -A UniqueAltReadCount \
    --native-pair-hmm-threads $NUM_PAIR_HMM_THREADS \
    --showHidden true
