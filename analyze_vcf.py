from glob import glob
from collections import Counter
from itertools import dropwhile
import subprocess
import csv

import numpy as np
import pandas as pd
import vcf


DATA_PATH = '/volume/cancer-mutect2/PCAWG_CLL_vcf_reindex'
file_names = glob(DATA_PATH + '/*.snv_mnv.vcf.gz')

white_list_donor_ids = set()

with open('./white_list_CLL_donor_ids.txt') as f:
    for line in f:
        donor_id = line.rstrip('\n')
        white_list_donor_ids.add(donor_id)

donor_file_name_map = {}
with open('./PCAWG_CLL_vcf_manifest.tsv') as f:
    csv_reader = csv.reader(f, delimiter='\t')
    header = next(csv_reader)
    donor_index = header.index('donor_id/donor_count')
    file_name_index = header.index('file_name')
    for row in csv_reader:
        if 'snv_mnv' not in row[file_name_index]:
            continue
        donor_id = row[donor_index]
        if donor_id not in white_list_donor_ids:
            continue
        if row[donor_index] in donor_file_name_map:
            donor_file_name_map[row[donor_index]].append(row[file_name_index])
        else:
            donor_file_name_map[row[donor_index]] = [ row[file_name_index]  ]

donor_ids = list(donor_file_name_map.keys())
# for debug purpose
# print(donor_file_name_map['DO6350'])
# for donor_id in donor_ids:
#     print(len(donor_file_name_map[donor_id]))

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
low_support = []
oxog_fail = []
bseq = []
bpcr = []
germ = []
germovlp = []
normal_panel = []
remapfail = []
for donor_id in donor_ids:
    consensus_file_name = ''
    dkfz_file_name = ''
    MUSE_file_name = ''
    svcp_file_name = ''
    broad_mutect_file_name = ''
    for file_name in donor_file_name_map[donor_id]:
        if 'consensus' in file_name:
            consensus_file_name = file_name
        elif 'dkfz' in file_name:
            dkfz_file_name = file_name
        elif 'MUSE' in file_name:
            MUSE_file_name = file_name
        elif 'svcp' in file_name:
            svcp_file_name = file_name
        elif 'mutect' in file_name:
            broad_mutect_file_name = file_name
    variant_counter = Counter()
    for vcf_file_name in [dkfz_file_name, MUSE_file_name, svcp_file_name, broad_mutect_file_name]:
        vcf_reader = vcf.Reader(open(DATA_PATH + '/' + vcf_file_name, 'rb'))
        for record in vcf_reader:
            for alt in record.ALT:
                # print(record.POS)
                variant_counter.update([(record.CHROM, record.POS, record.REF, str(alt))])
    total_varaints_of_four_pipelines.append(len(variant_counter))
    for key, count in dropwhile(lambda key_count: key_count[1] >= 2, variant_counter.most_common()):
        del variant_counter[key]
    # compare with current consensus
    consensus_set = set()
    consensus_vcf_reader = vcf.Reader(open(DATA_PATH + '/' + consensus_file_name, 'rb'))
    filter_counts = Counter()


    # print('number of SNVs that are low support: {}'.format(filter_counts['LOWSUPPORT']))
    # print('number of SNVs with failed OxOGfilter: {}'.format(filter_counts['OXOGFAIL']))
    # print('number of SNVs with seuqencing bias: {}'.format(filter_counts['bSeq']))
    # print('number of SNVs with PCR bias: {}'.format(filter_counts['bPcr']))
    # print('number of SNVs that are with insufficient somatic evidence in 1000 genome variants: {}'.format(filter_counts['GERM1000G']))
    # print('number of SNVs of Overlaps germline haplotype call: {}'.format(filter_counts['GERMOVLP']))
    # print('number of SNVs presented in panel of normals: {}'.format(filter_counts['NORMALPANEL']))
    # print('number of SNVs of no longer seen under remapping: {}'.format(filter_counts['REMAPFAIL']))


    
    for record in consensus_vcf_reader:
        for alt in record.ALT:
            consensus_set.add((record.CHROM, record.POS, record.REF, str(alt)))
            if record.FILTER is not None:
                for filter in record.FILTER:
                    filter_counts.update([filter])
    consensus_num_of_SNPs.append(len(consensus_set))
    all_qualified_variants = set(variant_counter)
    expected_2_plus_num_of_SNPs.append(len(all_qualified_variants))

    low_support.append(filter_counts['LOWSUPPORT'])
    oxog_fail.append(filter_counts['OXOGFAIL'])
    bseq.append(filter_counts['bSeq'])
    bpcr.append(filter_counts['bPcr'])
    germ.append(filter_counts['GERM1000G'])
    germovlp.append(filter_counts['GERMOVLP'])
    normal_panel.append(filter_counts['NORMALPANEL'])
    remapfail.append(filter_counts['REMAPFAIL'])


    # result = subprocess.run(['/volume/cancer-mutect2/bcftools-1.10.2/bcftools', 'stats', DATA_PATH + '/' + consensus_file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # for line in result.stdout.decode('utf-8').split('\n'):
    #    # print(line)
    #    if 'SN\t0\tnumber of SNPs:' in line:
    #        consensus_num_of_SNPs.append((int(line.split(':')[1].lstrip('\t'))))


    result = subprocess.run(['/volume/cancer-mutect2/bcftools-1.10.2/bcftools', 'stats', DATA_PATH + '/' + consensus_file_name, DATA_PATH + '/' + dkfz_file_name ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # for debug purpose
    # print(result.stderr)
    # print(result.stdout)
    for line in result.stdout.decode('utf-8').split('\n'):
        # print(line)
        if 'SN\t1\tnumber of SNPs:' in line:
            dkfz_num_of_SNPs.append(int(line.split(':')[1].lstrip('\t')))
        elif 'SN\t2\tnumber of SNPs:' in line:
            consensus_dkfz_common_num_of_SNPs.append(int(line.split(':')[1].lstrip('\t')))
            break
   
    result = subprocess.run(['/volume/cancer-mutect2/bcftools-1.10.2/bcftools', 'stats', DATA_PATH + '/' + consensus_file_name, DATA_PATH + '/' + MUSE_file_name ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in result.stdout.decode('utf-8').split('\n'):
        if 'Number of records' in line:
            print(line)
        if 'SN\t1\tnumber of SNPs:' in line:
            MUSE_num_of_SNPs.append(int(line.split(':')[1].lstrip('\t')))
        elif 'SN\t2\tnumber of SNPs:' in line:
            consensus_MUSE_common_num_of_SNPs.append(int(line.split(':')[1].lstrip('\t')))
            break

    result = subprocess.run(['/volume/cancer-mutect2/bcftools-1.10.2/bcftools', 'stats', DATA_PATH + '/' + consensus_file_name, DATA_PATH + '/' + svcp_file_name ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in result.stdout.decode('utf-8').split('\n'):
        if 'SN\t1\tnumber of SNPs:' in line:
            svcp_num_of_SNPs.append(int(line.split(':')[1].lstrip('\t')))
        elif 'SN\t2\tnumber of SNPs:' in line:
            consensus_svcp_common_num_of_SNPs.append(int(line.split(':')[1].lstrip('\t')))
            break

    result = subprocess.run(['/volume/cancer-mutect2/bcftools-1.10.2/bcftools', 'stats', DATA_PATH + '/' + consensus_file_name, DATA_PATH + '/' + broad_mutect_file_name ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in result.stdout.decode('utf-8').split('\n'):
        if 'SN\t1\tnumber of SNPs:' in line:
            broad_mutect_num_of_SNPs.append(int(line.split(':')[1].lstrip('\t')))
        elif 'SN\t2\tnumber of SNPs:' in line:
            consensus_broad_mutect_common_num_of_SNPs.append(int(line.split(':')[1].lstrip('\t')))
            break

print(len(total_varaints_of_four_pipelines))
print('Average total number of SNPs in four pipelines: {}'.format(np.average(total_varaints_of_four_pipelines)))
print('Average number of SNPs expected in the consensus set: {}'.format(np.average(expected_2_plus_num_of_SNPs)))
print('Average number of SNPs in consensus set: {}'.format(np.average(consensus_num_of_SNPs)))
print('Average number of SNPs called by dkfz but not in consensus set before filtering: {}'.format(np.average(dkfz_num_of_SNPs)))
print('Average number of SNPs called by MUSE but not in consensus set before filtering: {}'.format(np.average(MUSE_num_of_SNPs)))
print('Average number of SNPs called by svcp but not in consensus set before filtering: {}'.format(np.average(svcp_num_of_SNPs)))
print('Average number of SNPs called by broad mutect but not in consensus set: {}'.format(np.average(broad_mutect_num_of_SNPs)))
print('Average number of SNPs called by dkfz and in the consensus set : {}'.format(np.average(consensus_dkfz_common_num_of_SNPs)))
print('Average number of SNPs called by MUSE and in the consensus set : {}'.format(np.average(consensus_MUSE_common_num_of_SNPs)))
print('Average number of SNPs called by svcp and in the consensus set : {}'.format(np.average(consensus_svcp_common_num_of_SNPs)))
print('Average number of SNPs called by broad mutect and in the consensus set: {}'.format(np.average(consensus_broad_mutect_common_num_of_SNPs)))

print('Average:')
print('LOWSUPPORT', np.average(low_support))
print('OXOGFAIL', np.average(oxog_fail))
print('bSeq', np.average(bseq))
print('bPcr', np.average(bpcr))
print('GERM1000G', np.average(germ))
print('GERMOVLP', np.average(germovlp))
print('NORMALPANEL', np.average(normal_panel))
print('REMAPFAIL', np.average(remapfail))

print('standard deviation')
print('LOWSUPPORT', np.std(low_support))
print('OXOGFAIL', np.std(oxog_fail))
print('bSeq', np.std(bseq))
print('bPcr', np.std(bpcr))
print('GERM1000G', np.std(germ))
print('GERMOVLP', np.std(germovlp))
print('NORMALPANEL', np.std(normal_panel))
print('REMAPFAIL', np.std(remapfail))



df = pd.DataFrame(data={
        'total number of variants of four pipelines': total_varaints_of_four_pipelines,
        'expected_2_plus_num_of_SNPs': expected_2_plus_num_of_SNPs,
        'Number of SNPs in consensus set': consensus_num_of_SNPs,
        'Number of SNPs called by dkfz but not in consensus set': dkfz_num_of_SNPs,
        'Number of SNPs called by MUSE but not in consensus set': MUSE_num_of_SNPs,
        'Number of SNPs called by svcp but not in consensus set': svcp_num_of_SNPs,
        'Number of SNPs called by broad mutect but not in consensus set': broad_mutect_num_of_SNPs,
        'Number of SNPs called by dkfz and in the consensus set': consensus_dkfz_common_num_of_SNPs,
        'Number of SNPs called by MUSE and in the consensus set': consensus_MUSE_common_num_of_SNPs,
        'Number of SNPs called by svcp and in the consensus set': consensus_svcp_common_num_of_SNPs,
        'Number of SNPs called by broad mutect and in the consensus set': consensus_broad_mutect_common_num_of_SNPs
    })

df.to_pickle('./result.pkl')
