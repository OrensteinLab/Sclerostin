import pandas as pd
import Bio
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from Bio.Seq import Seq
from helper_functions import *
import re
from itertools import zip_longest

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip_longest(s1, s2))

def IUPAC(seq1, seq2):
    new_seq = ''
    if len(seq1) == len(seq2):
        for i in range(len(seq1)):
            if seq1[i] == seq2[i]:
                new_seq += seq1[i]
            elif (seq1[i] == 'A' and seq2[i] == 'G') or (seq2[i] == 'A' and seq1[i] == 'G'):
                new_seq += 'R'
            elif (seq1[i] == 'C' and seq2[i] == 'T') or (seq2[i] == 'C' and seq1[i] == 'T'):
                new_seq += 'Y'
            elif (seq1[i] == 'G' and seq2[i] == 'C') or (seq2[i] == 'G' and seq1[i] == 'C'):
                new_seq += 'S'
            elif (seq1[i] == 'A' and seq2[i] == 'T') or (seq2[i] == 'A' and seq1[i] == 'T'):
                new_seq += 'W'
            elif (seq1[i] == 'G' and seq2[i] == 'T') or (seq2[i] == 'G' and seq1[i] == 'T'):
                new_seq += 'K'
            elif (seq1[i] == 'A' and seq2[i] == 'C') or (seq2[i] == 'A' and seq1[i] == 'C'):
                new_seq += 'M'
            else:
                new_seq += 'N'
    return new_seq

def quality_merge(seq1, seq2, quality_seq1, quality_seq2):
    new_seq = ''
    if len(seq1) == len(seq2):
        for i in range(len(seq1)):
            if (seq1[i] == seq2[i]) or (quality_seq1[i] >= quality_seq2[i]):
                new_seq += seq1[i]
            else:
                new_seq += seq2[i]
    return new_seq

def merge_fw_rv(fw, rv, fw_quality, rv_quality):#, min_alignment=10, max_alignment=30):
    mismatch = []
    for i in range(10, 30):
        mismatch.append(hamming_distance(fw[len(fw) - i:], rv[:i]) / i)
    merge_pos = mismatch.index(min(mismatch)) + 10
    #iupac_merge = IUPAC(fw[len(fw) - merge_pos:], rv[:merge_pos])
    q_merge = quality_merge(fw[len(fw) - merge_pos:], rv[:merge_pos], fw_quality[len(fw) - merge_pos:], rv_quality[:merge_pos])
    merge_seq = fw[:len(fw) - merge_pos] + q_merge + rv[merge_pos:]
    return merge_seq


def align_seq(seq, WTnuc, i, j, k, rv=False):
    alignments = pairwise2.align.globalms(WTnuc[:len(seq)], seq, i, j, k, k)
    relevant_seq = alignments[0][1]
    if rv:
        relevant_seq = complementary_DNA(relevant_seq)
        aa_seq = translate_nuc_to_aa(relevant_seq[len(relevant_seq) % 3:])
    else:
        aa_seq = translate_nuc_to_aa(relevant_seq)
    stop_codon = aa_seq.find('*')
    unknown_codon = aa_seq.find('?')
    if stop_codon == -1 and unknown_codon == -1:
        return aa_seq


def slice_seq(seq, isReverse=0):
    start = 0
    relevant_seq = ''
    if isReverse:
        start_motif = re.findall(r'(?=([A,C,G,T]TACG|G[A,C,G,T]ACG|GT[A,C,G,T]CG|GTA[A,C,G,T]GC|GTAC[A,C,G,T]C|GTACG[A,C,G,T]))', seq)
    else:
        start_motif = re.findall(r'(?=([A,C,G,T]AAGG|C[A,C,G,T]AGG|CA[A,C,G,T]GG|CAA[A,C,G,T]GG|CAAG[A,C,G,T]G|CAAGG[A,C,G,T]))', seq)
    if start_motif:
        position = seq.find(start_motif[0])
        if position != -1:
            relevant_seq = seq[position:]
            if isReverse:
                relevant_seq = complementary_DNA(relevant_seq)
                start = len(relevant_seq) % 3
    return relevant_seq[start:]

def nuc_to_aa(seq):
    dna = Seq(seq)
    aa_seq = dna.translate()
    #aa_seq = translate_nuc_to_aa(seq)
    stop_codon = aa_seq.find('*')
    unknown_codon = aa_seq.find('X')
    if stop_codon != -1:
        return 'stop codon'
    elif unknown_codon != -1:
        return 'unknown codon'
    else:  # If the seq is valid
        return aa_seq
    #else:
        #if isReverse:
            #align_seq(seq, WTnuc, i=2, j=-6, k=-4, rv=1)
        #else:
            #align_seq(seq, WTnuc, i=2, j=-10, k=-6, rv=0)


def sort_seq(seq, WTnuc, isReverse=0):
    start = 0
    if isReverse:
        start_motif = re.findall(r'(?=([A,C,G,T]TACG|G[A,C,G,T]ACG|GT[A,C,G,T]CG|GTA[A,C,G,T]GC|GTAC[A,C,G,T]C|GTACG[A,C,G,T]))', seq)
    else:
        start_motif = re.findall(r'(?=([A,C,G,T]AAGG|C[A,C,G,T]AGG|CA[A,C,G,T]GG|CAA[A,C,G,T]GG|CAAG[A,C,G,T]G|CAAGG[A,C,G,T]))', seq)
    if start_motif:
        position = seq.find(start_motif[0])
        if position != -1:
            relevant_seq = seq[position:]
            if isReverse:
                relevant_seq = complementary_DNA(relevant_seq)
                start = len(relevant_seq) % 3
            aa_seq = translate_nuc_to_aa(relevant_seq[start:])
            stop_codon = aa_seq.find('*')
            unknown_codon = aa_seq.find('?')
            if stop_codon == -1 and unknown_codon == -1:  # If the seq is valid
               return aa_seq
            else:
                if isReverse:
                    align_seq(seq, WTnuc, i=2, j=-6, k=-4, rv=1)
                else:
                    align_seq(seq, WTnuc, i=2, j=-10, k=-6, rv=0)


def sort_data(seq_dict,isReverse=0):
    start = 0
    start_motif = []
    valid = pd.Series()
    not_valid = pd.Series()
    not_seq = pd.Series()
    for name, seq in seq_dict.items():
        if isReverse:
            start_motif = re.findall(
                r'(?=([A,C,G,T]TACG|G[A,C,G,T]ACG|GT[A,C,G,T]CG|GTA[A,C,G,T]GC|GTAC[A,C,G,T]C|GTACG[A,C,G,T]))', seq)
        else:
            start_motif = re.findall(
                r'(?=([A,C,G,T]AAGG|C[A,C,G,T]AGG|CA[A,C,G,T]GG|CAA[A,C,G,T]GG|CAAG[A,C,G,T]G|CAAGG[A,C,G,T]))', seq)
        if not start_motif:
            not_seq = not_seq.append(pd.Series(seq), ignore_index=True)
        else:
            position = seq.find(start_motif[0])
            if position == -1:
                not_seq = not_seq.append(pd.Series(seq), ignore_index=True)
            else:
                relevant_seq = seq[position:]
                if isReverse == 1:
                    relevant_seq = complementary_DNA(relevant_seq)
                    start = len(relevant_seq) % 3
                aa_seq = translate_nuc_to_aa(relevant_seq[start:])
                stop_codon = aa_seq.find('*')
                unknown_codon = aa_seq.find('?')
                if stop_codon == -1 and unknown_codon == -1:    # If the seq is valid
                    valid = valid.append(pd.Series(seq[position:]), ignore_index=True)
                else:
                    not_valid = not_valid.append(pd.Series(seq[position:]), ignore_index=True)
    return [not_seq, not_valid, valid]


def sort_to_mut(aa_seq, WTaa):
    if len(aa_seq) == len(WTaa):
        count_diff, diff, mut_pos = compare_seq_toWT(aa_seq, WTaa)
        return aa_seq, count_diff, diff, mut_pos


def sort_seqs_to_mut(seq_dict, start_motif, nuc_length, WTnuc, WTaa, isReverse=0, onlyseq=0):
    '''
    The function receives dictionary of sequences and returns dictionary sorted by the  valid sequences into mutations
    :param seq_dict: dictionary containing all the sequences
    :param start_motif: string- after the motif the desired sequence will begin
    :param nuc_length: int- the length of the desired sequence
    :param WTaa: string- The amino acid sequence of the WT
    :return: dictionary of sequences- the key: The name of the seq
                                      the values: aa_seq - amino acid seq, count_diff- number of differences,
                                       diff- string , mut_pos- list of mutation positions
    '''
    start = 0
    not_valid = 0
    ali_no_valid = 0
    not_exist = 0
    mutation_dict = dict()
    for name, seq in seq_dict.items():
        position = seq.find(start_motif)
        if position == -1:
            not_exist += 1
            print(seq)
        if position != -1: #and position + nuc_length <= len(seq):
            relevant_seq = seq[position:]
            if isReverse == 1:
                start = len(relevant_seq) - (int(len(relevant_seq) / 3) * 3) #start code from here
                relevant_seq = complementary_DNA(relevant_seq)
            aa_seq = translate_nuc_to_aa(relevant_seq[start:])
            stop_codon = aa_seq.find('*')
            unknown_codon = aa_seq.find('?')
            if stop_codon == -1 and unknown_codon == -1:    # If the seq is valid
                if isReverse == 1:
                    count_diff, diff, mut_pos = compare_seq_toWT(aa_seq[::-1], WTaa[::-1])
                else:
                    count_diff, diff, mut_pos = compare_seq_toWT(aa_seq, WTaa)
                if onlyseq == 1:
                    mutation_dict[name] = aa_seq
                else:
                    mutation_dict[name] = aa_seq, count_diff, diff, mut_pos
            else:
                alignments = pairwise2.align.globalms(WTnuc[:len(seq[position:])], seq[position:], 8, -10, -6, -4.8)
                relevant_seq = complementary_DNA(alignments[0][1])
                aa_seq = translate_nuc_to_aa(relevant_seq[len(relevant_seq) % 3:])
                #aa_seq = translate_nuc_to_aa(relevant_seq)
                #aa_seq = translate_nuc_to_aa(relevant_seq[start:])
                stop_codon = aa_seq.find('*')
                unknown_codon = aa_seq.find('?')
                if stop_codon == -1 and unknown_codon == -1:  # If the seq is valid
                    if isReverse == 1:
                        count_diff, diff, mut_pos = compare_seq_toWT(aa_seq[::-1], WTaa[::-1])
                    else:
                        count_diff, diff, mut_pos = compare_seq_toWT(aa_seq, WTaa)
                    if onlyseq == 1:
                        mutation_dict[name] = aa_seq
                    else:
                        mutation_dict[name] = aa_seq, count_diff, diff, mut_pos
                else:
                    ali_no_valid += 1
                    #print(seq)
                    #print(aa_seq)

                not_valid += 1
    print(not_valid)
    print(ali_no_valid)
    return mutation_dict


def sort_mut_by_number(df, num_of_mut, col_name):
    '''
    The function receives a dictionary of sequences, a number of mutations and a desired name
    and returns a series of read-counts for the mutations
    :param mut_dict: dictionary of sequences
    :param num_of_mut: int - number of mutation
    :param col_name: string - the desired name of the column
    :return: a series of read-counts for the mutations with the number of mutations
    '''
    #df = pd.DataFrame.from_records(data, columns=['seq', 'Number of Mutations', 'mutations', col_name])
    if num_of_mut < 5:
        ind_relevant_mut = df['Number of Mutations'] == num_of_mut
    else:
        ind_relevant_mut = df['Number of Mutations'] >= num_of_mut
    relevant_mut = df.loc[ind_relevant_mut, col_name].value_counts()
    return relevant_mut


