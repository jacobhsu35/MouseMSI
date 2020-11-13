#!/bin/bash

DIR=/project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/
INPUT_FILE_DIR=/project/GP1/u3710062/AI_SHARE/rawdata/Mouse_MSI/

mkdir -p ${DIR}/mut_ratio

echo 'Pair VCF FilterVCF NormVCF NormVCF_SNP NormVCF_INDEL PASS_VCF' > ${DIR}/mut_ratio/TN_pair.txt

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
	TN_DIR="${DIR}/${Tumor}_${i}"
	VCF="${TN_DIR}/${Tumor}_${i}_Mutect2.vcf.gz"
	FilterVCF="${TN_DIR}/${Tumor}_${i}_Mutect2.filtered.vcf.gz"
	NormVCF="${TN_DIR}/${Tumor}_${i}_Mutect2.norm.vcf"			# Normalized VCF
	NormVCF_SNP="${TN_DIR}//${Tumor}_${i}_Mutect2.norm.snp.vcf.gz"	# Normalized VCF (SNP)
	NormVCF_INDEL="${TN_DIR}/${Tumor}_${i}_Mutect2.norm.indel.vcf.gz"	# Normalized VCF (INDEL)
	PASS_VCF="${TN_DIR}/${Tumor}_${i}_Mutect2.norm.PASS.vcf.gz"		# Normalized VCF (only PASS or ".")
	### Count Variants
	VCF_N=`zgrep "BaseQRankSum=" $VCF | wc | awk '{print $1}'`
	FilterVCF_N=`zgrep "BaseQRankSum=" $FilterVCF | wc | awk '{print $1}'`
	NormVCF_N=`zgrep "BaseQRankSum=" $NormVCF | wc | awk '{print $1}'`
	NormVCF_SNP_N=`zgrep "BaseQRankSum=" $NormVCF_SNP | wc | awk '{print $1}'`
	NormVCF_INDEL_N=`zgrep "BaseQRankSum=" $NormVCF_INDEL | wc | awk '{print $1}'`
	PASS_VCF_N=`zgrep "BaseQRankSum=" $PASS_VCF | wc | awk '{print $1}'`
	Count="${Tumor}_${i} ${VCF_N} ${FilterVCF_N} ${NormVCF_N} ${NormVCF_SNP_N} ${NormVCF_INDEL_N} ${PASS_VCF_N}"
	echo "${Count}" >> ${DIR}/mut_ratio/${TODAY}_TN_pair.txt
	done
 done</project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/MouseMSI/INPUT/Mouse_MSI_TN_table.txt
