# Pipeline

## Prerequisite

- Python >= 3.5

## Usage

- `run_mutect2.sh`
  - An example is available [here](https://github.com/Jacob-s-Lab/CML/blob/12fc6b5071c5f747f1b90fe295d90410bc3f47c6/somatic_run.sh).

## Download Data

- Data Source:
  - Google Cloud Bucket of GATK Best Practice (somatic-b37): https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37
  - FTP sever of GATK

``` shell
# reference genome sequence
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.dict
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta.fai

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz
# the downloaded vcf file for dbSNP need to be decompressed first
gzip -d dbsnp_138.hg19.vcf.gz
gzip -d dbsnp_138.hg19.vcf.idx.gz
```

## Reference Links
- [MuSE](https://github.com/danielfan/MuSE)
