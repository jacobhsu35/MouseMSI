from glob import glob
import subprocess
import csv
import numpy as np
import pandas as pd

DATA_PATH = '/volume/cancer-mutect2/PCAWG_CLL_vcf_reindex'
file_names = glob(DATA_PATH + '/*.snv_mnv.vcf.gz')

donor_file_name_map = {}
with open('./PCAWG_CLL_vcf_manifest.tsv') as f:
    csv_reader = csv.reader(f, delimiter='\t')
    header = next(csv_reader)
    donor_index = header.index('donor_id/donor_count')
    file_name_index = header.index('file_name')
    for row in csv_reader:
        if 'snv_mnv'not in row[file_name_index]:
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

    result = subprocess.run(['/volume/cancer-mutect2/bcftools-1.10.2/bcftools', 'stats', DATA_PATH + '/' + consensus_file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in result.stdout.decode('utf-8').split('\n'):
        # print(line)
        if 'SN\t0\tnumber of SNPs:' in line:
            consensus_num_of_SNPs.append((int(line.split(':')[1].lstrip('\t'))))


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

print('Average number of SNPs in consensus set: {}'.format(np.average(consensus_num_of_SNPs)))
print('Average number of SNPs called by dkfz but not in consensus set: {}'.format(np.average(dkfz_num_of_SNPs)))
print('Average number of SNPs called by MUSE but not in consensus set: {}'.format(np.average(MUSE_num_of_SNPs)))
print('Average number of SNPs called by svcp but not in consensus set: {}'.format(np.average(svcp_num_of_SNPs)))
print('Average number of SNPs called by broad mutect but not in consensus set: {}'.format(np.average(broad_mutect_num_of_SNPs)))
print('Average number of SNPs called by dkfz and in the consensus set: {}'.format(np.average(consensus_dkfz_common_num_of_SNPs)))
print('Average number of SNPs called by MUSE and in the consensus set: {}'.format(np.average(consensus_MUSE_common_num_of_SNPs)))
print('Average number of SNPs called by svcp and in the consensus set: {}'.format(np.average(consensus_svcp_common_num_of_SNPs)))
print('Average number of SNPs called by broad mutect and in the consensus set: {}'.format(np.average(consensus_broad_mutect_common_num_of_SNPs)))

df = pd.DataFrame(data={
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
