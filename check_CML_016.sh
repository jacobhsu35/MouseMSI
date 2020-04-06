VCF_GZ_FILE=/work1/hsiaoyi0504/result/CML_hg19/TW-CML-M-016/Mutect2.vcf.gz
BCF_TOOLS=/pkg/biology/BCFtools/BCFtools_v1.10.2/bin/bcftools

# check location of variant
zcat $VCF_GZ_FILE | grep -P "chr2\t220098039"
zcat $VCF_GZ_FILE | grep -P "chr4\t71394929"
zcat $VCF_GZ_FILE | grep -P "chr8\t109241332"
zcat $VCF_GZ_FILE | grep -P "chr9\t134501557"
zcat $VCF_GZ_FILE | grep -P "chr17\t7324318"
zcat $VCF_GZ_FILE | grep -P "chr21\t36231791"
zcat $VCF_GZ_FILE | grep -P "chr22\t50355105"
zcat $VCF_GZ_FILE | grep -P "chrX\t39923044"
zcat $VCF_GZ_FILE | grep -P "chrX\t39931673"
zcat $VCF_GZ_FILE | grep -P "chrX\t82763798"

# check location of gene
$BCF_TOOLS view --no-header -r chr2:220094479-220101391 $VCF_GZ_FILE # ANKZF1 gene, https://www.genecards.org/cgi-bin/carddisp.pl?gene=ANKZF1
$BCF_TOOLS view --no-header -r chr4:71384257-71398459 $VCF_GZ_FILE # AMTN gene, https://www.genecards.org/cgi-bin/carddisp.pl?gene=AMTN
$BCF_TOOLS view --no-header -r chr8:109213445-109447562 $VCF_GZ_FILE # EIF3E gene, https://www.genecards.org/cgi-bin/carddisp.pl?gene=EIF3E
$BCF_TOOLS view --no-header -r chr9:134452157-134615461 $VCF_GZ_FILE # RAPGEF1,
$BCF_TOOLS view --no-header -r chr17:7323679-7324951 $VCF_GZ_FILE # SPEM1,
$BCF_TOOLS view --no-header -r chr21:36160098-37376965 $VCF_GZ_FILE # RUNX1,
$BCF_TOOLS view --no-header -r chr22:50354143-50357728  $VCF_GZ_FILE # PIM3,
$BCF_TOOLS view --no-header -r chrX:39909068-40036582  $VCF_GZ_FILE # BCOR,
$BCF_TOOLS view --no-header -r chrX:82763269-82764775 $VCF_GZ_FILE # POU3F4,
