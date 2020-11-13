# Pipeline

## Prerequisite

- Python >= 3.5

## Usage

### Predefine INPUT 
- `sample list and locations.xlsx`
  - Sample info from PI

- `Mouse_MSI_TN_table.txt`
  - This table is to define the T/N pairs as the input


MouseMSI_location_whole_exons.bed is from:
- `/Users/jacob_imac/Nextcloud/TAIWAN/NTU_research/DISEASES/Mouse_MSI/From_PI/sample\ list\ and\ locations.xlsx`

- `cl7_cMNR-panel_regions_cons.bed`
  - capture regions (BAIT)
  - /Users/jacob_imac/Nextcloud/TAIWAN/NTU_research/DISEASES/Mouse_MSI/From_PI/Ref_盈瑄論文/capture\ probe/cMNR\ probe\ design/Roche\ design/cl7_cMNR-panel_all_results/cl7_cMNR-panel_regions_cons.bed

- `cl7_cMNR-panel_capture_targets.bed`
  - capture regions (Target)
  - /Users/jacob_imac/Nextcloud/TAIWAN/NTU_research/DISEASES/Mouse_MSI/From_PI/Ref_盈瑄論文/capture\ probe/cMNR\ probe\ design/Roche\ design/cl7_cMNR-panel_all_results/Selection_Results/cl7_cMNR-panel_capture_targets.bed

### Bash @NCHC
- `Mouse_MSI.sh `
  - This bash is to run all pairs with GATK pipeline and VEP annotation.

- `run_mutect2.sh` 
  - This bash is to call mutations from FASTQ to VCF based on Mutect2 & HaplotypeCaller for each sample, respectively. 
  - An example is available [here](https://github.com/Jacob-s-Lab/CML/blob/12fc6b5071c5f747f1b90fe295d90410bc3f47c6/somatic_run.sh).

- `run_vep.sh`
  - To run VEP @ NCHC
  - Both VCF and Tab output format

- `run_vep_online.sh`
  - Query and to run VEP online

- `Calculate_MSI_change.sh`
  - To count the number of mutations
    - VCF=           Mutect2.vcf.gz			# Mutect2 default output
    - FilterVCF=     Mutect2.filtered.vcf.gz		# Filtered VCF (--min-allele-fraction 0.01 ; --unique-alt-read-count 20)
    - NormVCF=       Mutect2.norm.vcf			# Normalized VCF
    - NormVCF_SNP=   Mutect2.norm.snp.vcf.gz	   	# Normalized VCF (SNP)
    - NormVCF_INDEL= Mutect2.norm.indel.vcf.gz	 	# Normalized VCF (INDEL)
    - PASS_VCF=      Mutect2.norm.PASS.vcf.gz	  	# Normalized VCF (only PASS or ".")
  - OUTPUT=`TN_pair.txt`
