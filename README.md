# CLL

## Prerequisite
- Python >= 3.5

## Download Data

- Data Source:
  - Google Cloud Bucket of GATK Best Practice (somatic-b37): https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-b37
  - FTP sever of GATK

``` shell
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.dict
wget https://storage.googleapis.com/gatk-best-practices/somatic-b37/Homo_sapiens_assembly19.fasta.fai

wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz
```

## Reference Links
- [MuSE](https://github.com/danielfan/MuSE)
