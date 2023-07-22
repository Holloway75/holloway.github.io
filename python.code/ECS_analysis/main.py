import copy
import os
import numpy as np
from pyckmeans import CKmeans
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
import data_prepare
import fst_calculation
import plot
from addresss import *
import pandas as pd
from PoissonMixture import PoissonMixture
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.cluster import KMeans
import multiprocessing
from matplotlib_venn import venn3


def cor_matrix(input_df_, pc_num):
    input_df = copy.deepcopy(input_df_)
    sample_index_list = input_df.index.tolist()
    a = np.array(input_df[['PC%d' % i for i in range(1, pc_num+1)]])
    b = np.ones((a.shape[1], a.shape[0]))
    result = np.sqrt((np.dot(a**2, b) + np.dot(a**2, b).T - 2 * np.dot(a, a.T)))
    df_cor_matrx = pd.DataFrame(columns=sample_index_list, index=sample_index_list, dtype='float64', data=result)
    df_cor_matrx.fillna(0, inplace=True)
    return df_cor_matrx


def gene_freq_by_subset(input_df_, df_cor_matrix):
    input_df = copy.deepcopy(input_df_)
    sample_index_list = input_df.index.tolist()
    df_gene_filtered_result = pd.DataFrame(index=sample_index_list, columns=Auto_list, dtype='int8')
    for arr in sample_index_list:
        nearest_1000_list = df_cor_matrix.loc[arr].sort_values().index.tolist()[0:1000]
        df_tmp = input_df[[(i in nearest_1000_list) for i in input_df.index]]
        df_tmp = df_tmp[~df_tmp['gene'].isna()]
        gene_list = []
        for sample in df_tmp.index:
            gene_list += df_tmp.loc[sample, 'gene'].split(":")
        df_gene_filtered_result.loc[arr] = [gene_list.count(i) for i in Auto_list]
    return df_gene_filtered_result


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


def poisson_mixture_analysis(data, gene):
    columns = ['n_components', 'parameters', 'weights', 'bic']
    df = pd.DataFrame(columns=columns)
    for i in np.arange(1,8):
        mx = PoissonMixture(n_components=i).fit(data)
        df.loc[i, 'n_components'] = mx.n_components
        df.loc[i, 'parameters'] = ":".join(mx.params.reshape(-1).astype('str'))
        df.loc[i, 'weights'] = ":".join(mx.weights.reshape(-1).astype('str'))
        df.loc[i, 'bic'] = mx.bic(data)
    if gene == 'HBA1/HBA2':
        df.to_csv('D:\我的坚果云\ECS_1.6w_samples\gene_distribution_pc1\HBA1_HBA2.poisson_mixture.csv', encoding='utf_8_sig')
    else:
        df.to_csv('D:\我的坚果云\ECS_1.6w_samples\gene_distribution_pc1\%s.poisson_mixture.csv' % gene, encoding='utf_8_sig')


def poisson(lambda_, k, weight):
    lambda_, k, weight = float(lambda_), float(k), float(weight)
    log_prob = -lambda_ + k * np.log(lambda_) - sum([np.log(i) for i in np.arange(1, k+1)])
    return weight * 16090 * np.exp(log_prob)


def get_sample_name(dir):
    os.chdir(dir)
    fam_list = os.listdir()
    df = pd.DataFrame(columns=['fid', 'id_1', 'name_1', 'id_2', 'name_2'])
    for i in fam_list:
        os.chdir(dir + '\\' + i)
        sam_list = os.listdir()
        df.loc[i, 'fid'] = sam_list[0].split("-")[0]
        df.loc[i, 'id_1'] = sam_list[0].split("-")[1]
        df.loc[i, 'name_1'] = sam_list[0].split("-")[2]
        if len(sam_list) == 1:
            pass
        elif len(sam_list) == 2:
            df.loc[i, 'id_2'] = sam_list[1].split("-")[1]
            df.loc[i, 'name_2'] = sam_list[1].split("-")[2]
        else:
            print(i)
            raise ValueError
    return df


def plot_gene2(input_df, ax, ns, cut_line=1/200):
    pre_df = data_prepare.data2plot_gene(input_df, cut_line)
    if 'G6PD' in pre_df.columns:
        pre_df.drop(labels='G6PD', axis=1, inplace=True)
    y = pre_df.loc['total']
    gene_num = len(y)
    x = np.arange(gene_num)

    ax.barh(x, y, height=0.7, color='#1f77b4', align='center', tick_label=pre_df.columns.tolist())
    ax.set_xscale("log")
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    ax.axvline(1 / 200, color='r', label='Carrier frequency=1/200')
    ax.set_title('%d filtered autosomal genes in %s area' % (gene_num, ns))
    ax.legend()
    ax.set_xlabel('Carrier frequency')
    ax.set_ylabel('Gene')


if __name__ == '__main__':
    os.chdir(r'D:\我的坚果云\ECS_1.6w_samples')
    fig,ax = plt.subplots()
    df = pd.read_excel('ven3.xlsx')
    a = df['博圣'].dropna().tolist()
    b = df['协和'].dropna().tolist()
    c = df['课题'].dropna().tolist()
    venn3(subsets=[set(a), set(b), set(c)], set_labels=('博圣', '协和', '课题'))
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    plt.show()



    # df = pd.read_csv('gene_distribution_pc1/genes.different_cf.pc1.csv', index_col=0)
    # for i in df.index:
    #     cf_list = df.loc[i, 'cf_list'].split(":")
    #     df.loc[i, 'highest_cf'] = max([float(i) for i in cf_list])
    #     df.loc[i, 'w of subgroup with highest cf'] = df.loc[i, 'weights'].split(":")[cf_list.index(str(df.loc[i, 'highest_cf']))]
    # df.to_csv('gene_distribution_pc1/genes.different_cf.pc1.csv', index=True, encoding='utf_8_sig')
