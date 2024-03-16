import math
import pandas as pd
import numpy as np
import os
from sort_sclerostin import *
import time
import multiprocessing as mp
from functools import partial
startTime = time.time()

# Sclerostin library NGS- variables
WTaa = 'QGWQAFKNDATEIIPELGEYPEPPPELENNKTMNRAENGGRPPHHPFETKDVSEYSCRELHFTRYVTDGPCRSAKPVTELVCSGQCGPARLLPNAIGRGKWWRPSGPDF' \
       'RCIPDRYRAQRVQLLCPGGEAPRARKVRLVASCKCKRLTRFHNQSELKDFGTEAARPQKGRKPRPRARSAKANQAELENAY'
WTnucFW = 'CAAGGGTGGCAAGCGTTTAAAAATGACGCTACCGAAATCATCCCAGAGCTTGGTGAGTACCCCGAGCCCCCTCCAGAGCTAGAAAACAATAAAACAATGAATAGAGCT' \
        'GAGAATGGCGGCCGTCCCCCTCATCACCCCTTTGAAACTAAGGATGTCTCTGAATACTCCTGTAGGGAATTACACTTCACCAGGTACGTCACAGACGGCCCGTGTAGG' \
        'TCAGCAAAGCCGGTAACCGAGCTAGTATGTTCTGGTCAGTGTGGTCCTGCACGTTTGCTGCCAAATGCGATTGGTCGTGGAAAATGGTGGAGGCCATCCGGCCCGGAC' \
        'TTCCGTTGCATACCGGACCGTTACCGTGCACAACGTGTTCAGTTGCTGTGCCCAGGAGGGGAGGCCCCGAGGGCAAGGAAAGTTAGGCTTGTAGCATCTTGTAAGTGT' \
        'AAGAGGTTGACACGTTTTCATAACCAGAGTGAACTAAAAGACTTTGGTACAGAAGCAGCTAGACCTCAAAAGGGCCGTAAACCACGTCCAAGGGCCAGAAGTGCTAAA' \
        'GCTAACCAAGCCGAATTAGAGAACGCGTAC'
WTnucRV = complementary_DNA(WTnucFW)
start_code_fw = 'CAAGGGT'
start_code_rv = 'GTACGC'
nuc_length = len(WTnucFW)

mut_dict = {1: pd.DataFrame(), 2: pd.DataFrame(), 3: pd.DataFrame()}
mut_numbers_dict = {0: pd.DataFrame(), 1: pd.DataFrame(), 2: pd.DataFrame(), 3: pd.DataFrame(), 4: pd.DataFrame(), 5: pd.DataFrame()} #5- 5 or bigger


for name in mut_dict.keys():
    print(name)
    seq_dict_fw, q_dict_fw = fastq_reader('../' + 'LIB' + str(name) + 'FW' + '.fastq')
    seq_dict_rv, q_dict_rv = fastq_reader('../' + 'LIB' + str(name) + 'RV' + '.fastq')
    mut_dict[name] = pd.DataFrame({'fw': pd.Series(seq_dict_fw), 'rv': pd.Series(seq_dict_rv), 'fw_quality': pd.Series(q_dict_fw), 'rv_quality': pd.Series(q_dict_rv)})
    print(len(mut_dict[name]))
    mut_dict[name] = mut_dict[name].dropna()
    pool = mp.Pool(processes=6)
    df = pool.map(partial(slice_seq, isReverse=0), mut_dict[name]['fw'])
    mut_dict[name]['fw_seq'] = df
    df = pool.map(partial(slice_seq, isReverse=1), mut_dict[name]['rv'])
    mut_dict[name]['rv_seq'] = df

    mut_dict[name] = mut_dict[name].dropna()
    print(len(mut_dict[name]))
    all_pairs = []
    for pair in zip(mut_dict[name]['fw_seq'], mut_dict[name]['rv_seq'], mut_dict[name]['fw_quality'], mut_dict[name]['rv_quality']):
        all_pairs.append(pair)
    merged = pool.starmap(merge_fw_rv, all_pairs)
    mut_dict[name]['merged_seq'] = merged
    df = pool.map(nuc_to_aa, mut_dict[name]['merged_seq'])
    mut_dict[name]['aa_seq'] = df
    mut_dict[name].to_csv('all data merged lib aa' + str(name) + '.csv')
    pool.close()
    print('stop codon')
    print(len(mut_dict[name][mut_dict[name]['aa_seq'] == 'stop codon']))
    print('unknown codon')
    print(len(mut_dict[name][mut_dict[name]['aa_seq'] == 'unknown codon']))
    print('valid')
    print(len(mut_dict[name]) - len(mut_dict[name][mut_dict[name]['aa_seq'] == 'stop codon']) - len(mut_dict[name][mut_dict[name]['aa_seq'] == 'unknown codon']))

total = []

for name in mut_dict.keys():
    all_data = pd.read_csv('all data merged lib aa' + str(name) + '.csv')
    all_data = all_data.dropna()
    series = all_data['aa_seq'].apply(sort_to_mut, WTaa=WTaa)
    series = series.dropna()
    total.append(len(series))
    columns = ['seq', 'Number of Mutations', 'mutations' + str(name), 'positions' + str(name)]
    df = pd.DataFrame([[a, b, c, d] for a, b, c, d in series.values], columns=columns)
    for num_of_mut in mut_numbers_dict.keys():
        mut_numbers_dict[num_of_mut] = pd.concat([mut_numbers_dict[num_of_mut],
                                              sort_mut_by_number(df, num_of_mut, 'mutations' + str(name))], axis=1, sort=True)
        print('lib' + str(name) + '#mutations' + str(num_of_mut) + ':' + str(mut_numbers_dict[num_of_mut]['mutations' + str(name)].sum()))
print(total)


WT1 = int(mut_numbers_dict[0]['mutations1'])
WT2 = int(mut_numbers_dict[0]['mutations2'])
WT3 = int(mut_numbers_dict[0]['mutations3'])

for i in mut_numbers_dict.keys():
    mut_numbers_dict[i] = mut_numbers_dict[i].dropna()
    mut_numbers_dict[i]['freq_mutations1'] = mut_numbers_dict[i]['mutations1'] / total[0]
    mut_numbers_dict[i]['freq_mutations2'] = mut_numbers_dict[i]['mutations2'] / total[1]
    mut_numbers_dict[i]['freq_mutations3'] = mut_numbers_dict[i]['mutations3'] / total[2]

    mut_numbers_dict[i]['freq_WT_mutations1'] = mut_numbers_dict[i]['mutations1'] / WT1
    mut_numbers_dict[i]['freq_WT_mutations2'] = mut_numbers_dict[i]['mutations2'] / WT2
    mut_numbers_dict[i]['freq_WT_mutations3'] = mut_numbers_dict[i]['mutations3'] / WT3

    mut_numbers_dict[i]['2/1_WT'] = np.log2(mut_numbers_dict[i]['freq_WT_mutations2'] / mut_numbers_dict[i]['freq_WT_mutations1'])
    mut_numbers_dict[i]['3/1_WT'] = np.log2(mut_numbers_dict[i]['freq_WT_mutations3'] / mut_numbers_dict[i]['freq_WT_mutations1'])
    mut_numbers_dict[i]['3/2_WT'] = np.log2(mut_numbers_dict[i]['freq_WT_mutations3'] / mut_numbers_dict[i]['freq_WT_mutations2'])

    df = pd.concat([mut_numbers_dict[i]['2/1_WT'], mut_numbers_dict[i]['3/1_WT'], mut_numbers_dict[i]['3/2_WT'],
                      mut_numbers_dict[i]['mutations1'], mut_numbers_dict[i]['mutations2'], mut_numbers_dict[i]['mutations3'],
                      mut_numbers_dict[i]['freq_mutations1'], mut_numbers_dict[i]['freq_mutations2'],
                      mut_numbers_dict[i]['freq_mutations3']], axis=1, sort=True)
    df.to_csv('all data by variants ' + str(i) + ' mutations' + '.csv')
    test = build_all_possible_mut_df(WTaa)
    test = pd.concat([test, df], axis=1, sort=True)
    test = test.sort_values(by=['Pos', 'Mutation'])
    test.to_csv('all data by variants to heat map ' + str(i) + ' mutations' + '.csv')
