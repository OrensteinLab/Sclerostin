import math
import pandas as pd
import numpy as np
import os
from sort_sclerostin import *
import matplotlib.pyplot as plt

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

#paper heatmaps
import seaborn as sns

significant_variants_lib21 = pd.read_excel(open('data/significant variants.xlsx', 'rb'), sheet_name='ER2_1')
significant_variants_lib21 = significant_variants_lib21[significant_variants_lib21['p-adj'] < 0.05]
significant_variants_lib32 = pd.read_excel(open('data/significant variants.xlsx', 'rb'), sheet_name='ER3_2')
significant_variants_lib32 = significant_variants_lib32[significant_variants_lib32['p-adj'] < 0.05]

test = build_all_possible_mut_df(WTaa)
test = pd.concat([test, significant_variants_lib21['2/1_WT'], significant_variants_lib32['3/2_WT']], axis=1, sort=True)
test = test.sort_values(by=['Pos', 'Mutation'])
df = test.iloc[:, :]


def heatmap_df(df, title='', x_title='', y_title='', color='coolwarm', cbar_title='', figure_size=(7, 7), save_path=None):
    '''
    The function creates a bar plot showing val_counts with the specified parameters.
    :param df: A pandas data frame with the values to be plotted
    :param title: string, default is empty string. Title of the figure to be plotted
    :param x_title: string, default is empty string. Title (label) of the x axis.
    :param y_title: string, default is empty string. Title (label) of the y axis.
    :param color: string, default is 'coolwarm'. The color of the bars to be plotted
    :param cbar_title: string, default is empty string. Title (label) of the heat map scale.
    :param figure_size: tuple with 2 numbers, default is (7,7). The size of the figure.
    :param save_path: string, default is None. The full path of the file to which the figure should be saved.
           if save_path is None, the figure will be showed to the screen.

    The function creates heat map showing df with the specified parameters
    '''
    plt.figure(figsize=figure_size)
    ax = sns.heatmap(df, cmap=color, cbar_kws={'label': cbar_title}, vmin=-6, vmax=6)
    cbar = ax.collections[0].colorbar

    cbar.ax.yaxis.label.set_fontweight('bold')
    cbar.ax.tick_params(axis='y', labelsize=16, width=2)  # Adjust label size and width as needed
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontweight('bold')
        tick.set_family('Arial')  # Set font family to Arial

    # Set font properties
    cbar.set_label(cbar_title, fontsize=22, fontname='Arial')
    plt.grid()
    #plt.title(title, fontsize=24, fontweight='bold')
    plt.xlabel(x_title, fontsize=22, fontweight='bold', fontname='Arial')
    plt.ylabel(y_title, fontsize=22, fontweight='bold', fontname='Arial')
    #plt.xlim(0, len(df))
    #plt.yticks(np.arange(len(range(190))), [x+1 for x in range(191)], fontsize=8, rotation=0, fontweight='bold')

    num_ticks = len(df) // 5 + 1
    yticks = [1] + [i * 5 for i in range(1, num_ticks)]
    ytick_positions = np.linspace(0, len(df), len(yticks))
    plt.yticks(ytick_positions, yticks, fontsize=14, rotation=0, fontweight='bold', fontname='Arial')
    plt.gca().set_yticks(np.arange(len(df)), minor=True)
    plt.gca().yaxis.grid(True, which='minor', color='gray', linestyle='-', linewidth=1)


    plt.xticks(np.arange(len(range(20))),
               ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'], fontsize=14, fontweight='bold', rotation=0, fontname='Arial')
    ax.invert_yaxis()
    ax.axhline(y=0, color='k', linewidth=2)
    ax.axhline(y=df.shape[0], color='k', linewidth=2)
    ax.axvline(x=0, color='k', linewidth=2)
    ax.axvline(x=df.shape[1], color='k', linewidth=2)
    ax.tick_params(bottom=False)
    if save_path:
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
    else:
        plt.show()
    plt.close('all')

data = np.asarray(df['1/2_WT']).reshape(190, 20)

data_df = pd.DataFrame(data)
#df_all_transposed = data_df.transpose()

heatmap_df(data_df, title='', x_title='Amino acid', y_title='Position', figure_size=(6, 20),
           cbar_title='log2 enrichment ratio', save_path= 'data/log2_left' + '.jpg')

data = np.asarray(df['3/2_WT']).reshape(190, 20)

data_df = pd.DataFrame(data)
#df_all_transposed = data_df.transpose()

heatmap_df(data_df, title='', x_title='Amino acid', y_title='Position', figure_size=(6, 20),
           cbar_title='log2 enrichment ratio', save_path= 'data/log2_middle' + '.jpg')

test = build_all_possible_mut_df(WTaa)
merged = pd.concat([ significant_variants_lib21['2/1_WT'], significant_variants_lib32['3/2_WT']], axis=1, sort=True)
merged = merged.dropna()
test = pd.concat([test, merged], axis=1, sort=True)
test = test.sort_values(by=['Pos', 'Mutation'])
df = test.iloc[:, :]

data = np.asarray(df['1/2_WT']).reshape(190, 20)

data_df = pd.DataFrame(data)
#df_all_transposed = data_df.transpose()

heatmap_df(data_df, title='', x_title='Amino acid', y_title='Position', figure_size=(6, 20),
           cbar_title='log2 enrichment ratio', save_path= 'data/log2_right' + '.jpg')
