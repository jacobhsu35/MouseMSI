OUTPUT_DIR=/work1/hsiaoyi0504/result/CML_hg19/TW-CML-M-018
VCF_GZ_FILE=$OUTPUT_DIR/Mutect2.vcf.gz
BCF_TOOLS=/pkg/biology/BCFtools/BCFtools_v1.10.2/bin/bcftools
ANNOVAR_PATH=/pkg/biology/ANNOVAR/ANNOVAR_20191024

$ANNOVAR_PATH/convert2annovar.pl -format vcf4old $VCF_GZ_FILE > $OUTPUT_DIR/Mutect2.annovar 2> /dev/null

# check location of variant

cat $OUTPUT_DIR/Mutect2.annovar | grep -P "3\t12632420"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "3\t47037005"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "3\t49896871"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "3\t52255972"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "5\t94852706"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "5\t180219340"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "7\t4050672"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "9\t96070669"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "10\t124392747"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "13\t110437498"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "15\t96877834"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "17\t7578265"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "17\t20109120"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "19\t16513262"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "19\t46511501"
cat $OUTPUT_DIR/Mutect2.annovar | grep -P "20\t18462374"

