#!/bin/bash

### This is to count the number of variant in different settings. 

TODAY=`date +%Y%m%d`
DIR=/project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/
INPUT_FILE_DIR=/project/GP1/u3710062/AI_SHARE/rawdata/Mouse_MSI/
OUTPUT=${DIR}/mut_ratio/${TODAY}_TN_pair.txt

mkdir -p ${DIR}/mut_ratio

echo 'Pair M_VCF M_VCF_F M_VCF_F_Norm M_VCF_F_Norm_SNV M_VCF_F_Norm_INDEL M_VCF_F_Norm_PASS H_NVCF H_NVCF_F H_NVCF_F_Norm H_NVCF_F_Norm_SNV H_NVCF_F_Norm_INDEL H_NVCF_F_Norm_PASS H_TVCF H_TVCF_F H_TVCF_F_Norm H_TVCF_F_Norm_SNV H_TVCF_F_Norm_INDEL H_TVCF_F_Norm_PASS ' > ${OUTPUT}

while read -r Tumor Ctrl_1 Ctrl_2 Ctrl_3 Ctrl_4 Ctrl_5 Ctrl_6 Ctrl_7;
 do
 array=($Ctrl_1 $Ctrl_2 $Ctrl_3 $Ctrl_4 $Ctrl_5 $Ctrl_6 $Ctrl_7)
 for i in "${array[@]}"
        do
	cd ${DIR}/mut_ratio
        wkdir=${DIR}/${Tumor}_${i}
        TODAY=`date +%Y%m%d`
        logfile=${DIR}/mut_ratio/${TODAY}_${Tumor}_${i}_mut_ratio.log
        exec 3<&1 4<&2
        exec >$logfile 2>&1
        set -euo pipefail
        set -x
	### Define input
	### Nomenclature --> Caller(M/H)_(VCF)_Filter(F)_(Norm)_(SNV/INDEL/PASS)_Count(C)SNV/INDEL
	### Only () will be used
	TN_DIR="${DIR}/${Tumor}_${i}"
	M_VCF="${TN_DIR}/${Tumor}_${i}_Mutect2.vcf.gz"
	M_VCF_F="${TN_DIR}/${Tumor}_${i}_Mutect2.filtered.vcf.gz"
	M_VCF_F_Norm="${TN_DIR}/${Tumor}_${i}_Mutect2.norm.vcf"			# Normalized M_VCF
	M_VCF_F_Norm_SNV="${TN_DIR}//${Tumor}_${i}_Mutect2.norm.snp.vcf.gz"		# Normalized M_VCF (SNV)
	M_VCF_F_Norm_INDEL="${TN_DIR}/${Tumor}_${i}_Mutect2.norm.indel.vcf.gz"	# Normalized M_VCF (INDEL)
	M_VCF_F_Norm_PASS="${TN_DIR}/${Tumor}_${i}_Mutect2.norm.PASS.vcf.gz"		# Normalized M_VCF (only PASS or ".")
	### VCF from HaplotypeCaller
	# Normal part
	H_NVCF="${TN_DIR}/${i}_marked.recal.pass1.haplotype.SnpIndel.vcf.gz"			#
	H_NVCF_F="${TN_DIR}/${i}_marked.recal.pass1.filtered.haplotype.SnpIndel.vcf"
	H_NVCF_F_Norm="${TN_DIR}/${i}_marked.recal.pass1.filtered.haplotype.SnpIndel.norm.vcf"
	H_NVCF_F_Norm_SNV="${TN_DIR}/${i}_marked.recal.pass1.filtered.haplotype.Snp.norm.vcf.gz"
	H_NVCF_F_Norm_INDEL="${TN_DIR}/${i}_marked.recal.pass1.filtered.haplotype.Indel.norm.vcf.gz"
	H_NVCF_F_Norm_PASS="${TN_DIR}/${i}_marked.recal.pass1.filtered.haplotype.SnpIndel.norm.PASS.vcf.gz"
	# Tumor part
	H_TVCF="${TN_DIR}/${Tumor}_marked.recal.pass1.haplotype.SnpIndel.vcf.gz"                    
        H_TVCF_F="${TN_DIR}/${Tumor}_marked.recal.pass1.filtered.haplotype.SnpIndel.vcf"
        H_TVCF_F_Norm="${TN_DIR}/${Tumor}_marked.recal.pass1.filtered.haplotype.SnpIndel.norm.vcf"
        H_TVCF_F_Norm_SNV="${TN_DIR}/${Tumor}_marked.recal.pass1.filtered.haplotype.Snp.norm.vcf.gz"
        H_TVCF_F_Norm_INDEL="${TN_DIR}/${Tumor}_marked.recal.pass1.filtered.haplotype.Indel.norm.vcf.gz"
        H_TVCF_F_Norm_PASS="${TN_DIR}/${Tumor}_marked.recal.pass1.filtered.haplotype.SnpIndel.norm.PASS.vcf.gz"
	### Count for T-N Variants
	M_VCF_C=`zgrep "BaseQRankSum=" $M_VCF | wc | awk '{print $1}'`
	M_VCF_F_C=`zgrep "BaseQRankSum=" $M_VCF_F | wc | awk '{print $1}'`
	M_VCF_F_Norm_C=`zgrep "BaseQRankSum=" $M_VCF_F_Norm | wc | awk '{print $1}'`
	M_VCF_F_Norm_SNV_C=`zgrep "BaseQRankSum=" $M_VCF_F_Norm_SNV | wc | awk '{print $1}'`
	M_VCF_F_Norm_INDEL_C=`zgrep "BaseQRankSum=" $M_VCF_F_Norm_INDEL | wc | awk '{print $1}'`
	M_VCF_F_Norm_PASS_C=`zgrep "BaseQRankSum=" $M_VCF_F_Norm_PASS | wc | awk '{print $1}'`
	### Count for HaplotypeCaller T and N
	H_NVCF_C=`zgrep "BaseQRankSum=" $H_NVCF | wc | awk '{print $1}'`
	H_NVCF_F_C=`zgrep "BaseQRankSum=" $H_NVCF_F | wc | awk '{print $1}'`
	H_NVCF_F_Norm_C=`zgrep "BaseQRankSum" $H_NVCF_F_Norm| wc | awk '{print $1}'`
	H_NVCF_F_Norm_SNV_C=`zgrep "BaseQRankSum" $H_NVCF_F_Norm_SNV| wc | awk '{print $1}'`
	H_NVCF_F_Norm_INDEL_C=`zgrep "BaseQRankSum" $H_NVCF_F_Norm_INDEL| wc | awk '{print $1}'`
	H_NVCF_F_Norm_PASS_C=`zgrep "BaseQRankSum" $H_NVCF_F_Norm_PASS| wc | awk '{print $1}'`
	H_TVCF_C=`zgrep "BaseQRankSum=" $H_TVCF | wc | awk '{print $1}'`
        H_TVCF_F_C=`zgrep "BaseQRankSum=" $H_TVCF_F | wc | awk '{print $1}'`
        H_TVCF_F_Norm_C=`zgrep "BaseQRankSum" $H_TVCF_F_Norm| wc | awk '{print $1}'`
        H_TVCF_F_Norm_SNV_C=`zgrep "BaseQRankSum" $H_TVCF_F_Norm_SNV| wc | awk '{print $1}'`
        H_TVCF_F_Norm_INDEL_C=`zgrep "BaseQRankSum" $H_TVCF_F_Norm_INDEL| wc | awk '{print $1}'`
        H_TVCF_F_Norm_PASS_C=`zgrep "BaseQRankSum" $H_TVCF_F_Norm_PASS| wc | awk '{print $1}'`
	###
	Count="${Tumor}_${i} ${M_VCF_C} ${M_VCF_F_C} ${M_VCF_F_Norm_C} ${M_VCF_F_Norm_SNV_C} ${M_VCF_F_Norm_INDEL_C} ${M_VCF_F_Norm_PASS_C} ${H_NVCF_C} ${H_NVCF_F_C} ${H_NVCF_F_Norm_C} ${H_NVCF_F_Norm_SNV_C} ${H_NVCF_F_Norm_INDEL_C} ${H_NVCF_F_Norm_PASS_C} ${H_TVCF_C} ${H_TVCF_F_C} ${H_TVCF_F_Norm_C} ${H_TVCF_F_Norm_SNV_C} ${H_TVCF_F_Norm_INDEL_C} ${H_TVCF_F_Norm_PASS_C}"
	echo "${Count}" >> ${OUTPUT}
	done
 done</project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/MouseMSI/INPUT/Mouse_MSI_TN_table.txt
### replace space to tab
sed -i 's/ /\t/g' ${OUTPUT}
### Description for each column 
echo "Pair: Tumor-Normal pair" >> ${OUTPUT}
echo "M_VCF: (Mutect2) Somatic mutation with default parameter " >> ${OUTPUT}
echo "M_VCF_F: (Mutect2) Somatic mutation with FilterMutectCalls function (min-allele-fraction=0.01, unique-alt-read-count=20)" >> ${OUTPUT}
echo "M_VCF_F_Norm:  M_VCF_F normalized by VT" >> ${OUTPUT}
echo "M_VCF_F_Norm_SNV: Only SNV of M_VCF_F_Norm" >> ${OUTPUT}
echo "M_VCF_F_Norm_INDEL: Only INDEL variant of M_VCF_F_Norm" >> ${OUTPUT}
echo "M_VCF_F_Norm_PASS: All PASS variant of M_VCF_F_Norm" >> ${OUTPUT}
echo "H_TVCF: (HaplotypeCaller) Germline mode on TUMOR part with default parameter" >> ${OUTPUT}
echo "H_TVCF_F: (HaplotypeCaller) Germline mode with VariantFiltration function (cluster-window-size=10;LowCoverage: DP < 10;VeryLowQual: QUAL < 30.0;LowQual: QUAL > 30.0 && QUAL < 50.0;LowQD: QD < 1.5) " >> ${OUTPUT}
echo "H_TVCF_F_Norm: H_TVCF_F normalized by VT" >> ${OUTPUT}
echo "H_TVCF_F_Norm_SNV: Only SNV of H_TVCF_F_Norm" >> ${OUTPUT}
echo "H_TVCF_F_Norm_INDEL: Only INDEL variant of H_TVCF_F_Norm" >> ${OUTPUT}
echo "H_TVCF_F_Norm_PASS: All PASS variant of H_TVCF_F_Norm" >> ${OUTPUT}
echo "H_NVCF: (HaplotypeCaller) Germline mode on NORMAL part with default parameter" >> ${OUTPUT}
echo "H_NVCF_F: (HaplotypeCaller) Germline mode with VariantFiltration function (cluster-window-size=10;LowCoverage: DP < 10;VeryLowQual: QUAL < 30.0;LowQual: QUAL > 30.0 && QUAL < 50.0;LowQD: QD < 1.5) " >> ${OUTPUT}
echo "H_NVCF_F_Norm: H_NVCF_F normalized by VT" >> ${OUTPUT}
echo "H_NVCF_F_Norm_SNV: Only SNV of H_NVCF_F_Norm" >> ${OUTPUT}
echo "H_NVCF_F_Norm_INDEL: Only INDEL variant of H_NVCF_F_Norm" >> ${OUTPUT}
echo "H_NVCF_F_Norm_PASS: All PASS variant of H_NVCF_F_Norm" >> ${OUTPUT}

