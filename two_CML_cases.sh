SCRIPT_PATH=/home/hsiaoyi0504/CLL
INPUT_FILE_DIR=/work1/hsiaoyi0504/data/CML

$SCRIPT_PATH/run_mutect2.sh $INPUT_FILE_DIR/TW-CML-M-016-BC_R1.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-016-BC_R2.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-016-D_R1.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-016-D_R2.fastq.bz2 \
    "@RG\tID:GP-TW-CML-M-016-BC\tSM:TW-CML-M-016-BC\tPL:illumina" \
    "@RG\tID:GP-TW-CML-M-016-D\tSM:TW-CML-M-016-D\tPL:illumina" \
    /work1/hsiaoyi0504/result/CML_hg19/TW-CML-M-016 \
    /work2/u00srx00/reference/GATK/human_g1k_v37/human_g1k_v37_decoy.fasta \
    /project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/2.8/b37/dbsnp_138.b37.vcf \
    /project/GP1/reference/Homo_sapiens/GATK/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz \
    /project/GP1/reference/Homo_sapiens/GATK/Mutect2/GetPileupSummaries/small_exac_common_3_b37.vcf.gz

$SCRIPT_PATH/run_mutect2.sh $INPUT_FILE_DIR/TW-CML-M-018-BC_R1.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-018-BC_R2.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-018-D_R1.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-018-D_R2.fastq.bz2 \
    "@RG\tID:GP-TW-CML-M-018-BC\tSM:TW-CML-M-018-BC\tPL:illumina" \
    "@RG\tID:GP-TW-CML-M-018-D\tSM:TW-CML-M-018-D\tPL:illumina" \
    /work1/hsiaoyi0504/result/CML_hg19/TW-CML-M-018 \
    /work2/u00srx00/reference/GATK/human_g1k_v37/human_g1k_v37_decoy.fasta \
    /project/GP1/u3710062/AI_SHARE/reference/GATK_bundle/2.8/b37/dbsnp_138.b37.vcf \
    /project/GP1/reference/Homo_sapiens/GATK/Mutect2/af-only-gnomad.raw.sites.b37.vcf.gz \
    /project/GP1/reference/Homo_sapiens/GATK/Mutect2/GetPileupSummaries/small_exac_common_3_b37.vcf.gz
