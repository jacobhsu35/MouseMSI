OUTPUT_DIR=/work1/hsiaoyi0504/result/CML_hg19/TW-CML-M-016
VCF_GZ_FILE=$OUTPUT_DIR/Mutect2.vcf.gz
BCF_TOOLS=/pkg/biology/BCFtools/BCFtools_v1.10.2/bin/bcftools
ANNOVAR_PATH=/pkg/biology/ANNOVAR/ANNOVAR_20191024

$ANNOVAR_PATH/convert2annovar.pl -format vcf4old $VCF_GZ_FILE > $OUTPUT_DIR/Mutect2.annovar 2> /dev/null 

# check location of variant
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "2\t220098039"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "4\t71394929"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "8\t109241332"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "9\t134501557"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "17\t7324318"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "21\t36231791"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "22\t50355105"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "X\t39923044"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "X\t39931673"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "X\t82763798"

