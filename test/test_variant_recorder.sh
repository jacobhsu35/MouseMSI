HTSLIB_PATH=/pkg/biology/HTSLIB/HTSLIB_v1.10.2/bin/
export PATH=$HTSLIB_PATH:$PATH

module load biology/Perl/default

/pkg/biology/Ensembl-VEP/Ensembl-VEP_v99.2/variant_recoder -i ./test_data/test2AMLvariants.vcf --grch37 --pretty

/pkg/biology/Ensembl-VEP/Ensembl-VEP_v99.2/variant_recoder -i ./test_data/test2AMLvariants_update_4.vcf --grch37 --pretty

/pkg/biology/Ensembl-VEP/Ensembl-VEP_v99.2/variant_recoder -i ./test_data/test2AMLvariants_update_4.vcf --grch37 --host asiadb.ensembl.org

/pkg/biology/Ensembl-VEP/Ensembl-VEP_v99.2/variant_recoder -i ./test_data/test2AMLvariants_update_4.vcf --host asiadb.ensembl.org
