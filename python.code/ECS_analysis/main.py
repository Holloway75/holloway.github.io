import copy
import os
os.environ["OMP_NUM_THREADS"] = '1'
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import data_prepare
import fst_calculation
import plot
from addresss import *
import pandas as pd
from PoissonMixture import PoissonMixture
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from pyckmeans import CKmeans


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
        df.to_csv('D:\我的坚果云\ECS_1.6w_samples\gene_distribution_pc1_new\HBA1_HBA2.poisson_mixture.csv', encoding='utf_8_sig')
    else:
        df.to_csv('D:\我的坚果云\ECS_1.6w_samples\gene_distribution_pc1_new\%s.poisson_mixture.csv' % gene, encoding='utf_8_sig')


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


if __name__ == '__main__':
    os.chdir(r'D:\我的坚果云\ECS_1.6w_samples')
    df = pd.read_excel('gene_freq_counts.xlsx', index_col=0)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0.15, 0.1, 0.75, 0.8])
    sns.lineplot(data=df, x='max_freq', y='y_%', ax=ax, color='black')
    ax.set_xscale('log')
    ax.set_xlim(0.0005, 0.04)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([1/1000, 1/200, 1/100, 1/50])
    ax.set_xticklabels(['1/1000', '1/200', '1/100', '1/50'])
    ax.set_yticks([0.25, 0.5, 0.6236, 0.7483, 0.9747])
    ax.set_yticklabels(['25%', '50%', '62.36%', '74.83%', '97.47%'])
    ax.vlines(x=0.01, ymin=0, ymax=0.6236, colors='blue', linestyles='dashed')
    ax.vlines(x=0.005, ymin=0, ymax=0.7483, colors='blue', linestyles='dashed')
    ax.vlines(x=0.001, ymin=0, ymax=0.9747, colors='blue', linestyles='dashed')
    ax.hlines(y=0.6236, xmin=0.01, xmax=0.04, linestyle='dashed', color='blue')
    ax.hlines(y=0.7483, xmin=0.005, xmax=0.04, linestyle='dashed', color='blue')
    ax.hlines(y=0.9747, xmin=0.001, xmax=0.04, linestyle='dashed', color='blue')
    ax.set_xlabel('Threshold of carrier frequency')
    ax.set_ylabel('Detection rate of variants')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.invert_xaxis()
    plt.show()

