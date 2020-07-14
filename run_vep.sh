VEP_CACHE_DIR=/pkg/biology/DATABASES/vep-j116831/Cache
VEP_PLUGIN_DIR=/pkg/biology/DATABASES/vep-j116831/Plugins
VEP_PLUGIN_DATA_DIR=/pkg/biology/DATABASES/vep-j116831/Data_for_plugins
VEP_PATH=/pkg/biology/Ensembl-VEP/Ensembl-VEP_v99.2
HTSLIB_PATH=/pkg/biology/HTSLIB/HTSLIB_v1.10.2/bin/
export PATH=$HTSLIB_PATH:$PATH

module load biology/Perl/default
# Don't need to run vt, since the test example is known to be fine

INPUT_VCF_PATH=$1
OUTPUT_VCF_PATH=$2

$VEP_PATH/vep --cache --offline \
    --cache_version 100 \
    --assembly GRCh37 \
    --port 3337 \
    --dir_plugins $VEP_PLUGIN_DIR \
    --dir_cache $VEP_CACHE_DIR \
    -i $INPUT_VCF_PATH \
    --vcf \
    -o $OUTPUT_VCF_PATH \
    --check_existing \
    --plugin LoFtool,$VEP_PLUGIN_DIR/LoFtool_scores.txt \
    --plugin ExACpLI,$VEP_PLUGIN_DIR/ExACpLI_values.txt \
    --plugin MPC,$VEP_PLUGIN_DATA_DIR/fordist_constraint_official_mpc_values_v2.txt.gz \
    --plugin LOVD \
    --plugin FlagLRG,$VEP_PLUGIN_DATA_DIR/list_LRGs_transcripts_xrefs.txt \
    --plugin FunMotifs,$VEP_PLUGIN_DATA_DIR/blood.funmotifs_sorted.bed.gz,fscore,dnase_seq \
    --plugin PostGAP,$VEP_PLUGIN_DATA_DIR/postgap_GRCh37.txt.gz,ALL \
    --plugin satMutMPRA,file=$VEP_PLUGIN_DATA_DIR/satMutMPRA_GRCh37_ALL.gz,cols=ALL \
    --fork 4

