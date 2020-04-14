from glob import glob
from collections import Counter
from itertools import dropwhile
import subprocess
import csv
import os

import numpy as np
import pandas as pd
import vcf


DATA_PATH = '/volume/cancer-mutect2/PCAWG_CLL_vcf_reindex'
PASS_ONLY_DATA_PATH = '/volume/cancer-mutect2/final_consensus_snv_indel_passonly/snv_mnv'
file_names = glob(DATA_PATH + '/*.snv_mnv.vcf.gz')


consensus_num_of_SNPs = []
dkfz_num_of_SNPs = []
MUSE_num_of_SNPs= []
svcp_num_of_SNPs = []
broad_mutect_num_of_SNPs = []
consensus_dkfz_common_num_of_SNPs = []
consensus_MUSE_common_num_of_SNPs = []
consensus_svcp_common_num_of_SNPs = []
consensus_broad_mutect_common_num_of_SNPs = []
expected_2_plus_num_of_SNPs = []
total_varaints_of_four_pipelines = []

for file_id in ['dc505248-ed04-4f77-a7c6-3fefbc5df27b', '08f7812b-0d74-42ba-985b-d0a027e8a80c', 
        '1193a9c4-5aab-4cd7-a690-60c96bd1172d', 'fa9996c3-b874-4424-a5a6-f1c7c0f42b9b']:

    print('file_id:', file_id)
    final_passonly_consensus_file_path = glob(PASS_ONLY_DATA_PATH + '/' + file_id + '*.snv_mnv.vcf.gz')[0]
    consensus_file_name = os.path.basename(glob(DATA_PATH + '/{}.consensus.*.somatic.snv_mnv.vcf.gz'.format(file_id))[0])
    dkfz_file_name = os.path.basename(glob(DATA_PATH + '/{}.dkfz-snvCalling*.somatic.snv_mnv.vcf.gz'.format(file_id))[0])
    MUSE_file_name = os.path.basename(glob(DATA_PATH + '/{}.MUSE_*.somatic.snv_mnv.vcf.gz'.format(file_id))[0])
    svcp_file_name = os.path.basename(glob(DATA_PATH + '/{}.svcp_*.somatic.snv_mnv.vcf.gz'.format(file_id))[0])
    broad_mutect_file_name = os.path.basename(glob(DATA_PATH + '/{}.broad-mutect*.somatic.snv_mnv.vcf.gz'.format(file_id))[0])
    num_variants_by_caller = Counter()
    variant_counter = Counter()
    for vcf_file_name in [dkfz_file_name, MUSE_file_name, svcp_file_name, broad_mutect_file_name]:
        vcf_reader = vcf.Reader(open(DATA_PATH + '/' + vcf_file_name, 'rb'))
        num = 0
        for record in vcf_reader:
            for alt in record.ALT:
                variant_counter.update([(record.CHROM, record.POS, record.REF, str(alt))])
                num_variants_by_caller.update([vcf_file_name])
                num += 1
        print('number of SNVs in {}: {}'.format(vcf_file_name, num))
    print('number of SNVs from four pipeline (union): {}'.format(len(variant_counter)))

    for key, count in dropwhile(lambda key_count: key_count[1] >= 2, variant_counter.most_common()):
        del variant_counter[key]
    print('number of SNVs from four pipeline (at least called by two different callers): {}'.format(len(variant_counter)))

    # compare with current consensus
    consensus_set = set()
    consensus_vcf_reader = vcf.Reader(open(DATA_PATH + '/' + consensus_file_name, 'rb'))
    filter_counts = Counter()
    num_callers_counts = Counter()
    caller_counts = Counter()
    caller_two_plus_counts = Counter()
    num_oxog_fail = 0
    num = 0
    for record in consensus_vcf_reader:
        if record.FILTER is not None:
            for filter in record.FILTER:
                filter_counts.update([filter])
        # if len(record.ALT) != 1:
        #    print('It should not happen')
        for alt in record.ALT:
            num += 1
            consensus_set.add((record.CHROM, record.POS, record.REF, str(alt)))
        num_callers = record.INFO['NumCallers']
        num_callers_counts.update([str(num_callers)])
        if num_callers != 1:
            for caller in record.INFO['Callers']:
                caller_counts.update([caller])
                caller_two_plus_counts.update([caller])
        else:
            caller_counts.update([record.INFO['Callers'][0]])

    print('------ Following lines were counted from {} file'.format(consensus_file_name))
    print('number of SNVs in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(num)) 
    print('number of SNVs called by dkfz pipeline in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(caller_counts['dkfz']))
    print('number of SNVs called by MUSE pipeline in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(caller_counts['muse']))
    print('number of SNVs called by svcp pipeline in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(caller_counts['sanger']))
    print('number of SNVs called by broad Mutect pipeline in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(caller_counts['broad']))
    print('number of SNVs called by dkfz pipeline will be kept after two plus voting: {}'.format(caller_two_plus_counts['dkfz']))
    print('number of SNVs called by MUSE pipeline will be kept after two plus voting: {}'.format(caller_two_plus_counts['muse']))
    print('number of SNVs called by svcp pipeline will be kept after two plus voting: {}'.format(caller_two_plus_counts['sanger']))
    print('number of SNVs called by broad Mutect pipeline will be kept after two plus voting: {}'.format(caller_two_plus_counts['broad']))
    # TODO
    # print('number of SNVs will be called')

    # for debug
    # print(filter_counts)

    ##FILTER=<ID=LOWSUPPORT,Description="Not called by enough callers in ensemble">
    ##FILTER=<ID=OXOGFAIL,Description="Failed OXOG oxidative artifact filter">
    ##FILTER=<ID=bSeq,Description="Sequencing Bias">
    ##FILTER=<ID=bPcr,Description="PCR Bias">
    ##FILTER=<ID=GERM1000G,Description="1000Genome variant with insufficient somatic evidence">
    ##FILTER=<ID=GERMOVLP,Description="Overlaps germline Haplotype call">
    ##FILTER=<ID=NORMALPANEL,Description="Presence in Panel of Normals">
    ##FILTER=<ID=REMAPFAIL,Description="Variant no longer seen under remapping">
    print('number of SNVs called by only one caller: {}'.format(num_callers_counts['1']))
    print('number of SNVs called by two callers: {}'.format(num_callers_counts['2']))
    print('number of SNVs called by three callers: {}'.format(num_callers_counts['3']))
    print('number of SNVs called by all four callers: {}'.format(num_callers_counts['4']))
    print('number of SNVs that are low support: {}'.format(filter_counts['LOWSUPPORT']))
    print('number of SNVs with failed OxOGfilter: {}'.format(filter_counts['OXOGFAIL']))
    print('number of SNVs with seuqencing bias: {}'.format(filter_counts['bSeq']))
    print('number of SNVs with PCR bias: {}'.format(filter_counts['bPcr']))
    print('number of SNVs that are with insufficient somatic evidence in 1000 genome variants: {}'.format(filter_counts['GERM1000G']))
    print('number of SNVs of Overlaps germline haplotype call: {}'.format(filter_counts['GERMOVLP']))
    print('number of SNVs presented in panel of normals: {}'.format(filter_counts['NORMALPANEL']))
    print('number of SNVs of no longer seen under remapping: {}'.format(filter_counts['REMAPFAIL']))
    
    disappered_variants = set(variant_counter) - consensus_set
    # print(disappered_variants)
    
    # from final pass only vcf file
    final_consensus_set = set()
    final_num_callers_counts = Counter()
    final_caller_counts = Counter()
    pass_only_consensus_vcf_reader = vcf.Reader(open(final_passonly_consensus_file_path, 'rb'))
    num = 0
    for record in pass_only_consensus_vcf_reader:
        for alt in record.ALT:
            num += 1
            final_consensus_set.add((record.CHROM, record.POS, record.REF, str(alt)))
        num_callers = record.INFO['NumCallers']
        final_num_callers_counts.update([str(num_callers)])
        for caller in record.INFO['Callers']:
            final_caller_counts.update([caller])

    print('------ Following lines were counted from final pass only {} file'.format(os.path.basename(final_passonly_consensus_file_path)))
    print('number of SNVs in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(num)) 
    print('number of SNVs called by dkfz pipeline in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(final_caller_counts['dkfz']))
    print('number of SNVs called by MUSE pipeline in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(final_caller_counts['muse']))
    print('number of SNVs called by svcp pipeline in *consensus.*.somatic.snv_mnv.vcf.gz file: {}'.format(final_caller_counts['sanger']))
    print('number of SNVs called by broad Mutect pipeline in *consensus.*.somatic.snv_mnv.vcf.gz file: {}\n'.format(final_caller_counts['broad']))
    # print('number of SNVs called by dkfz pipeline will be kept after two plus voting: {}'.format(final_caller_two_plus_counts['dkfz']))
    # print('number of SNVs called by MUSE pipeline will be kept after two plus voting: {}'.format(final_caller_two_plus_counts['muse']))
