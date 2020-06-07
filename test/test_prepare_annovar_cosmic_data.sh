ANNOVAR_PATH=/pkg/biology/ANNOVAR/ANNOVAR_20191127
DATA_PATH=/project/GP1/j116831/AI_Labs/data/COSMIC

./prepare_annovar_user.pl -dbtype cosmic $DATA_PATH/CosmicMutantExport_GRCh37_v91.tsv -vcf $DATA_PATH/CosmicCodingMuts_GRCh37_v91.normal.vcf > $DATA_PATH/cosmic_coding_GRCh37_v91.txt
./prepare_annovar_user.pl -dbtype cosmic $DATA_PATH/CosmicNCV_GRCh37_v91.tsv -vcf $DATA_PATH/CosmicNonCodingVariants_GRCh37_v91.normal.vcf > $DATA_PATH/cosmic_noncoding_GRCh37_v91.txt
