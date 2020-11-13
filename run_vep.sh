#!/bin/bash
#PBS -l select=1:ncpus=10
#PBS -q ntu192G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N VEP_SKLS518-thymus_SKLS518-kidney
#PBS -o /project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/SKLS518-thymus_SKLS518-kidney/20200812_SKLS518-thymus_SKLS518-kidney_VEP.out
#PBS -e /project/GP1/u3710062/AI_SHARE/GATK/Outputs/20200507_MouseMSI/SKLS518-thymus_SKLS518-kidney/20200812_SKLS518-thymus_SKLS518-kidney_VEP.err
#PBS -M jacobhsu@ntu.edu.tw
#PBS -m e
VEP_CACHE_DIR=/pkg/biology/DATABASES/vep-j116831/Cache
VEP_PLUGIN_DIR=/pkg/biology/DATABASES/vep-j116831/Plugins
VEP_PLUGIN_DATA_DIR=/pkg/biology/DATABASES/vep-j116831/Data_for_plugins
VEP_PATH=/pkg/biology/Ensembl-VEP/Ensembl-VEP_v99.2
HTSLIB_PATH=/pkg/biology/HTSLIB/HTSLIB_v1.10.2/bin/
export PATH=$HTSLIB_PATH:$PATH

module load biology/Perl/default
# Don't need to run vt, since the test example is known to be fine

INPUT_VCF_PATH=$OUTPUT_PATH/${TUMOR_ID}_${NORMAL_ID}_Mutect2.vcf.gz
OUTPUT_VCF_PATH=$OUTPUT_PATH/VEP_${TUMOR_ID}_${NORMAL_ID}

$VEP_PATH/vep --cache --offline \
    --cache_version 99 \
    --merged \
    --species mus_musculus \
    --assembly GRCm38 \
    --dir_plugins $VEP_PLUGIN_DIR \
    --dir_cache $VEP_CACHE_DIR \
    -i $INPUT_VCF_PATH \
    --vcf \
    -o ${OUTPUT_VCF_PATH}.vcf \
    --check_existing \
    --fork 4 \
    --force_overwrite 

$VEP_PATH/vep --cache --offline \
    --cache_version 99 \
    --merged \
    --species mus_musculus \
    --assembly GRCm38 \
    --dir_plugins $VEP_PLUGIN_DIR \
    --dir_cache $VEP_CACHE_DIR \
    -i $INPUT_VCF_PATH \
    --tab \
    -o ${OUTPUT_VCF_PATH}.txt \
    --check_existing \
    --fork 4 \
    --force_overwrite

#    --plugin LoFtool,$VEP_PLUGIN_DIR/LoFtool_scores.txt \

#    --plugin ExACpLI,$VEP_PLUGIN_DIR/ExACpLI_values.txt \
#    --plugin MPC,$VEP_PLUGIN_DATA_DIR/fordist_constraint_official_mpc_values_v2.txt.gz \
#    --plugin LOVD \
#    --plugin FlagLRG,$VEP_PLUGIN_DATA_DIR/list_LRGs_transcripts_xrefs.txt \
#    --plugin FunMotifs,$VEP_PLUGIN_DATA_DIR/blood.funmotifs_sorted.bed.gz,fscore,dnase_seq \
#    --plugin PostGAP,$VEP_PLUGIN_DATA_DIR/postgap_GRCh37.txt.gz,ALL \
#    --plugin satMutMPRA,file=$VEP_PLUGIN_DATA_DIR/satMutMPRA_GRCh37_ALL.gz,cols=ALL \
#    --fork 4

