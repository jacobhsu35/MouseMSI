#!/bin/bash
#PBS -l select=1:ncpus=10
#PBS -q ngs48G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N MouseMSI_TN
#PBS -o /project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/SKLS518_20200717_test/20200721_Mutect2.out
#PBS -e /project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/SKLS518_20200717_test/20200721_Mutect2.err
#PBS -M jacobhsu@ntu.edu.tw
#PBS -m e

SCRIPT_PATH=/project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/MouseMSI/
INPUT_FILE_DIR=/project/GP1/u3710062/AI_SHARE/rawdata/Mouse_MSI/
#INPUT_FILE_DIR=/project/GP1/u3710062/AI_SHARE/rawdata/Mouse_MSI/fastq_1st_samples/SKLS1345-kidney_S12_L001_R2_001.fastq.gz
#Ref=/project/GP1/u3710062/AI_SHARE/reference/GRCm38/GCA_000001635.8_GRCm38.p6_genomic.fna
Ref=/project/GP1/u3710062/AI_SHARE/reference/GRCm38_sanger/GRCm38_68.fa
DBSNP=/project/GP1/u3710062/AI_SHARE/reference/GRCm38_dbSNP/00-All.vcf.gz
BED=/project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/For_Picard/MouseMSI_location_whole_exons.interval_list

###/project/GP1/u3710062/AI_SHARE/rawdata/Mouse_MSI/*/SKLS1345-thymus*R1*fastq.gz
###/project/GP1/u3710062/AI_SHARE/rawdata/Mouse_MSI/*/SKLS1345-thymus*R2*fastq.gz


while read -r Tumor Ctrl_1 Ctrl_2 Ctrl_3 Ctrl_4 Ctrl_5 Ctrl_6 Ctrl_7;
 do
 array=($Ctrl_1 $Ctrl_2 $Ctrl_3 $Ctrl_4 $Ctrl_5 $Ctrl_6 $Ctrl_7)
 for i in "${array[@]}"
	do
	cd ${SCRIPT_PATH}
	mkdir -p ../${Tumor}_${i}
	wkdir=${SCRIPT_PATH}../${Tumor}_${i}
	cd ${wkdir}
	TODAY=`date +%Y%m%d`
	logfile=${wkdir}/${TODAY}_${Tumor}_${i}_run.log
	TODAY=`date +%Y%m%d`
	exec 3<&1 4<&2
	exec >$logfile 2>&1
	set -euo pipefail
	set -x   
	### To add PBS parameters and arguments
	echo '
#!/bin/bash
#PBS -l select=1:ncpus=10
#PBS -q ngs48G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N Mutect2_'${Tumor}_${i}'
#PBS -o '${wkdir}'/'${TODAY}'_'${Tumor}'_'${i}'_Mutect2.out
#PBS -e '${wkdir}'/'${TODAY}'_'${Tumor}'_'${i}'_Mutect2.err
#PBS -M jacobhsu@ntu.edu.tw
#PBS -m e
GATK_PATH=/pkg/biology/GATK/GATK_v4.1.8.0
PICARD_PATH=/pkg/biology/Picard/Picard_v2.18.11/picard.jar
BWA_PATH=/pkg/biology/BWA/BWA_v0.7.17
SAMTOOLS_PATH=/pkg/biology/SAMtools/SAMtools_v1.10/bin
BCFTOOLS_PATH=/pkg/biology/BCFtools/BCFtools_v1.10.2/bin
SVABA_PATH=/pkg/biology/SvABA/SvABA_v1.1.0/bin
TUMOR_FASTQ_1_PATH='${INPUT_FILE_DIR}'/*_samples/'${Tumor}'*R1*.fastq.gz
TUMOR_FASTQ_2_PATH='${INPUT_FILE_DIR}'/*_samples/'${Tumor}'*R2*.fastq.gz
NORMAL_FASTQ_1_PATH='${INPUT_FILE_DIR}'/*_samples/'${i}'*R1*.fastq.gz
NORMAL_FASTQ_2_PATH='${INPUT_FILE_DIR}'/*_samples/'${i}'*R2*.fastq.gz
TUMOR_ID='${Tumor}'
NORMAL_ID='${i}'
OUTPUT_PATH='${wkdir}'
REF_GENOME_PATH='${Ref}'
HUMAN_DBSNP_PATH='${DBSNP}'
INTERVAL='${BED}'
	' > ${wkdir}/run_${Tumor}_${Ctrl_1}_mutect2.sh
	### To add GATK preprocess and Mutect2 -L
	tail -n +37 ${SCRIPT_PATH}/run_mutect2.sh >> ${wkdir}/run_${Tumor}_${Ctrl_1}_mutect2.sh
	sed -i -e '190d' ${wkdir}/run_${Tumor}_${Ctrl_1}_mutect2.sh
	sed -i '189a\ -bamout $OUTPUT_PATH/Mutect2.bam > $OUTPUT_PATH/mutect2.log 2>&1' ${wkdir}/run_${Tumor}_${Ctrl_1}_mutect2.sh
	sed -i '189a\ -L ${INTERVAL}\ \\' ${wkdir}/run_${Tumor}_${Ctrl_1}_mutect2.sh
	chmod 755 ${wkdir}/run_${Tumor}_${Ctrl_1}_mutect2.sh
#	qsub ${wkdir}/run_${Tumor}_${Ctrl_1}_mutect2.sh
#	sleep 10s
	done
 done</project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/MouseMSI/Mouse_MSI_TN_table.txt
#$SCRIPT_PATH/run_mutect2.sh  \
#	$INPUT_FILE_DIR/*_sample/SKLS518-thymus_*_R1_001.fastq.gz \
#	$INPUT_FILE_DIR/SKLS518-thymus_*_R2_001.fastq.gz \
#	$INPUT_FILE_DIR/SKLS518-kidney_*_R1_001.fastq.gz \
#	$INPUT_FILE_DIR/SKLS518-kidney_*_R2_001.fastq.gz \
#	"SKLS518-thymus" \
#	"SKLS518-kidney" \
#	${SCRIPT_PATH}/../SKLS518_20200717_test \
#	${Ref} \
#	/project/GP1/u3710062/AI_SHARE/reference/GRCm38_dbSNP/00-All.vcf.gz \
#	/project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/For_Picard/MouseMSI_location_whole_exons.interval_list

#	/project/GP1/u3710062/AI_SHARE/reference/GRCm38_dbSNP/00-All.vcf.gz \
#	/project/GP1/u3710062/AI_SHARE/reference/GRCm38_dbSNP/00-All.vcf.gz 

#/project/GP1/u3710062/AI_SHARE/reference/GRCm38_sanger/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz 
#done</project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/MouseMSI/Mouse_MSI_TN_table.txt
