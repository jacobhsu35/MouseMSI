SCRIPT_PATH=/home/hsiaoyi0504/CLL
INPUT_FILE_DIR=/work1/hsiaoyi0504/data/CML

$SCRIPT_PATH/run_mutect2.sh $INPUT_FILE_DIR/TW-CML-M-016-BC_R1.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-016-BC_R2.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-016-D_R1.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-016-D_R2.fastq.bz2 \
    "@RG\tID:GP-TW-CML-M-016-BC\tSM:TW-CML-M-016-BC\tPL:illumina" \
    "@RG\tID:GP-TW-CML-M-016-D\tSM:TW-CML-M-016-D\tPL:illumina" \
    /work1/hsiaoyi0504/result/CML_hg19_2/TW-CML-M-016 \
    /home/hsiaoyi0504/ucsc.hg19.fasta \
    /work1/hsiaoyi0504/ref/hg19/dbsnp_138.hg19.vcf

$SCRIPT_PATH/run_mutect2.sh $INPUT_FILE_DIR/TW-CML-M-018-BC_R1.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-018-BC_R2.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-018-D_R1.fastq.bz2 \
    $INPUT_FILE_DIR/TW-CML-M-018-D_R2.fastq.bz2 \
    "@RG\tID:GP-TW-CML-M-018-BC\tSM:TW-CML-M-018-BC\tPL:illumina" \
    "@RG\tID:GP-TW-CML-M-018-D\tSM:TW-CML-M-018-D\tPL:illumina" \
    /work1/hsiaoyi0504/result/CML_hg19_2/TW-CML-M-018 \
    /home/hsiaoyi0504/ucsc.hg19.fasta \
    /work1/hsiaoyi0504/ref/hg19/dbsnp_138.hg19.vcf
