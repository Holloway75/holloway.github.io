import copy
from scipy.special import factorial
import numpy as np
from sklearn.cluster import KMeans
from addresss import *
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import poisson
from PoissonMixture import Poisson_Mixture


def cor_matrix(input_df_, pc_num):
    input_df = copy.deepcopy(input_df_)
    sample_index_list = input_df.index.tolist()
    a = np.array(input_df[['PC%d' % i for i in range(1, pc_num+1)]])
    b = np.ones((a.shape[1], a.shape[0]))
    result = np.sqrt((np.dot(a**2, b) + np.dot(a**2, b).T - 2 * np.dot(a, a.T)))
    df_cor_matrx = pd.DataFrame(columns=sample_index_list, index=sample_index_list, dtype='float64', data=result)
    df_cor_matrx.fillna(0, inplace=True)
    df_cor_matrx.to_csv('D:\cor_matirx\cor_matrix_pc%d.csv'%pc_num, index=True)


def gene_freq_by_subset(input_df_, pc_num):
    input_df = copy.deepcopy(input_df_)
    sample_index_list = input_df.index.tolist()
    df_cor_matrx = pd.read_csv('D:\cor_matirx\cor_matrix_pc%d.csv'%pc_num, index_col=0)
    df_gene_filtered_result = pd.DataFrame(index=sample_index_list, columns=Auto_list, dtype='int8')
    for arr in sample_index_list:
        nearest_1000_list = df_cor_matrx.loc[arr].sort_values().index.tolist()[0:1000]
        input_df.loc[arr, 'nearest_1000_list'] = ":".join(nearest_1000_list)
        df_tmp = input_df[[(i in nearest_1000_list) for i in input_df.index]]
        df_tmp = df_tmp[~df_tmp['gene'].isna()]
        gene_list = []
        for sample in df_tmp.index:
            gene_list += df_tmp.loc[sample, 'gene'].split(":")
        df_gene_filtered_result.loc[arr] = [gene_list.count(i) for i in Auto_list]
        print('pc%d, %f' % (pc_num, sample_index_list.index(arr)/len(sample_index_list)))
    input_df.to_csv('D:\cor_matirx\sample.nearest_1000_by-pc%d.csv' % pc_num, index=True, encoding='utf_8_sig')
    df_gene_filtered_result.to_csv('D:\cor_matirx\gene_filtered_result_pc%d.csv'% pc_num, index=True, encoding='utf_8_sig')


def max_gene_freq():
    index_list = ['not_grouped'] + ['pc%d'%i for i in range(1, 21)]
    col_list = Auto_list
    df = pd.DataFrame(dtype='float64', columns=col_list, index=index_list)

    df_sample = pd.read_csv('D:\我的坚果云\ECS_mid\sample.coord.csv', index_col='indivID')
    sample_amount = df_sample.shape[0]
    df_sample = df_sample[~df_sample['gene'].isna()]
    gene_list = []
    for i in df_sample.index:
        gene_list += df_sample.loc[i, 'gene'].split(":")
    for i in Auto_list:
        df.loc['not_grouped', i] = gene_list.count(i)/sample_amount

    for pc_num in range(1,21):
        df_nearest = pd.read_csv('sample.nearest_1000_by-pc%d.csv'%pc_num, index_col='indivID')
        keep_list = []
        set_list = []
        for i in df_nearest.index:
            set_tmp = set(df_nearest.loc[i, 'nearest_1000_list'].split(":"))
            if set_tmp not in set_list:
                keep_list.append(i)
                set_list.append(set_tmp)
        df_gene_freq = pd.read_csv('gene_filtered_result_pc%d.csv' % pc_num, index_col=0)
        df_gene_freq = df_gene_freq[[(i in keep_list) for i in df_gene_freq.index]]
        for gene in Auto_list:
            df.loc['pc%d'%pc_num, gene] = df_gene_freq[gene].max() / 1000

    # df = pd.read_csv('max_freq.csv', index_col=0)
    df = df.T
    df.to_csv('max_freq1.csv', index=True, encoding='utf_8_sig')


if __name__ == '__main__':
    df = pd.read_excel('D:/pigfile.xlsx
    print(df)