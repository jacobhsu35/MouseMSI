VCF_GZ_FILE=/work1/hsiaoyi0504/result/CML_hg19/TW-CML-M-018/Mutect2.vcf.gz
BCF_TOOLS=/pkg/biology/BCFtools/BCFtools_v1.10.2/bin/bcftools

zcat $VCF_GZ_FILE | grep -P "chr3\t12632420"
zcat $VCF_GZ_FILE | grep -P "chr3\t47037005"
zcat $VCF_GZ_FILE | grep -P "chr3\t49896871"
zcat $VCF_GZ_FILE | grep -P "chr3\t52255972"
zcat $VCF_GZ_FILE | grep -P "chr5\t94852706"
zcat $VCF_GZ_FILE | grep -P "chr5\t180219340"
zcat $VCF_GZ_FILE | grep -P "chr7\t4050672"
zcat $VCF_GZ_FILE | grep -P "chr9\t96070669"
zcat $VCF_GZ_FILE | grep -P "chr10\t124392747"
zcat $VCF_GZ_FILE | grep -P "chr13\t110437498"
zcat $VCF_GZ_FILE | grep -P "chr15\t96877834"
zcat $VCF_GZ_FILE | grep -P "chr17\t7578265"
zcat $VCF_GZ_FILE | grep -P "chr17\t20109120"
zcat $VCF_GZ_FILE | grep -P "chr19\t16513262"
zcat $VCF_GZ_FILE | grep -P "chr119\t46511501"
zcat $VCF_GZ_FILE | grep -P "chr119\t18462374"

# check location of gene
$BCF_TOOLS view --no-header -r chr3:12625100-12705725 $VCF_GZ_FILE # RAF1 gene, https://www.genecards.org/cgi-bin/carddisp.pl?gene=RAF1
$BCF_TOOLS view --no-header -r chr3:47021173-47051193 $VCF_GZ_FILE # NBEAL2 
$BCF_TOOLS view --no-header -r chr3:49895421-49907655 $VCF_GZ_FILE # CAMKV
$BCF_TOOLS view --no-header -r chr3:52255096-52273183 $VCF_GZ_FILE # TLR9
$BCF_TOOLS view --no-header -r chr5:94799599-94890711 $VCF_GZ_FILE # TTC37
$BCF_TOOLS view --no-header -r chr5:180217541-180242652 $VCF_GZ_FILE # MGAT1
$BCF_TOOLS view --no-header -r chr7:3341080-4308632 $VCF_GZ_FILE # SDK1
$BCF_TOOLS view --no-header -r chr9:95947198-96082854  $VCF_GZ_FILE # WNK2
$BCF_TOOLS view --no-header -r chr10:124320181-124403252 $VCF_GZ_FILE # DMBT1
$BCF_TOOLS view --no-header -r chr13:110406184-110438915 $VCF_GZ_FILE # IRS2
$BCF_TOOLS view --no-header -r chr15:96869157-96883492 $VCF_GZ_FILE # NR2F2
$BCF_TOOLS view --no-header -r chr17:7565097-7590863 $VCF_GZ_FILE # TP53
$BCF_TOOLS view --no-header -r chr17:19912657-20222339 $VCF_GZ_FILE # SPECC1
$BCF_TOOLS view --no-header -r chr19:16466050-16582896 $VCF_GZ_FILE # EPS15L1
$BCF_TOOLS view --no-header -r chr19:46498339-46524576 $VCF_GZ_FILE # CCDC61
$BCF_TOOLS view --no-header -r chr20:18447771-18465287 $VCF_GZ_FILE # POLR3F
