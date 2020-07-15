HTSLIB_PATH=/pkg/biology/HTSLIB/HTSLIB_v1.10.2/bin/
export PATH=$HTSLIB_PATH:$PATH

module load biology/Perl/default

/pkg/biology/Ensembl-VEP/Ensembl-VEP_v99.2/variant_recoder -i ./test_data/test2AMLvariants.vcf --grch37 --pretty > ./test_data_output/test2AMLvariants.variant_recoder.json
