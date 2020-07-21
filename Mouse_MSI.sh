#!/bin/bash
#PBS -l select=1:ncpus=10
#PBS -q nhri192G
#PBS -P MST108173
#PBS -W group_list=MST108173
#PBS -N SKLS518_20200717_test
#PBS -o /project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/SKLS518_20200717_test/20200721_Mutect2.out
#PBS -e /project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/SKLS518_20200717_test/20200721_Mutect2.err
#PBS -M jacobhsu@ntu.edu.tw
#PBS -m e

SCRIPT_PATH=/project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/MouseMSI/
INPUT_FILE_DIR=/project/GP1/u3710062/AI_SHARE/rawdata/Mouse_MSI/fastq_TEST_run/
#INPUT_FILE_DIR=/project/GP1/u3710062/AI_SHARE/rawdata/Mouse_MSI/fastq_1st_samples/SKLS1345-kidney_S12_L001_R2_001.fastq.gz
#Ref=/project/GP1/u3710062/AI_SHARE/reference/GRCm38/GCA_000001635.8_GRCm38.p6_genomic.fna
Ref=/project/GP1/u3710062/AI_SHARE/reference/GRCm38_sanger/GRCm38_68.fa

#sleep 120m &&

cd ${SCRIPT_PATH}
mkdir -p ../SKLS518_20200717_test
cd ../SKLS518_20200717_test

#$SCRIPT_PATH/mouse_mutect2.sh \
#$SCRIPT_PATH/SKLS518_20200717_test/20200716_run_mouse_mutect2.sh \
#$SCRIPT_PATH/SKLS518_20200717_test/20200719_run_mouse_mutect2.sh \
#$SCRIPT_PATH/SKLS518_20200717_test/20200720_run_mouse_mutect2.sh \
$SCRIPT_PATH/run_mutect2.sh  \
	$INPUT_FILE_DIR/SKLS518-thymus_*_R1_001.fastq.gz \
	$INPUT_FILE_DIR/SKLS518-thymus_*_R2_001.fastq.gz \
	$INPUT_FILE_DIR/SKLS518-kidney_*_R1_001.fastq.gz \
	$INPUT_FILE_DIR/SKLS518-kidney_*_R2_001.fastq.gz \
	"SKLS518-thymus" \
	"SKLS518-kidney" \
	${SCRIPT_PATH}/../SKLS518_20200717_test \
	${Ref} \
	/project/GP1/u3710062/AI_SHARE/reference/GRCm38_dbSNP/00-All.vcf.gz \
	/project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/For_Picard/MouseMSI_location_whole_exons.interval_list

#	/project/GP1/u3710062/AI_SHARE/reference/GRCm38_dbSNP/00-All.vcf.gz \
#	/project/GP1/u3710062/AI_SHARE/reference/GRCm38_dbSNP/00-All.vcf.gz 

#/project/GP1/u3710062/AI_SHARE/reference/GRCm38_sanger/mgp.v6.merged.norm.snp.indels.sfiltered.vcf.gz 

