import copy
import os
import math
import multiprocessing
import itertools
import matplotlib.ticker as ticker
import numpy as np
from scipy import linalg
from sklearn.cluster import KMeans
from addresss import *
import pandas as pd
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot as plt
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


def poisson_analysis(data, gene, pc_num):
    columns = ['n_components', 'parameters', 'weights', 'bic']
    df = pd.DataFrame(columns=columns)
    for i in np.arange(1,8):
        mx = Poisson_Mixture(n_components=i).fit(data)
        df.loc[i, 'n_components'] = mx.n_components
        df.loc[i, 'parameters'] = ":".join(mx.params.reshape(-1).astype('str'))
        df.loc[i, 'weights'] = ":".join(mx.weights.reshape(-1).astype('str'))
        df.loc[i, 'bic'] = mx.bic(data)
    if gene == 'HBA1/HBA2':
        df.to_csv('D:\我的坚果云\ECS_mid\gene_distribution_pc%d\HBA1_HBA2.poisson_mixture.csv'%pc_num, encoding='utf_8_sig')
    else:
        df.to_csv('D:\我的坚果云\ECS_mid\gene_distribution_pc%d\%s.poisson_mixture.csv' % (pc_num, gene), encoding='utf_8_sig')


def poisson(lambda_, k, weight):
    return weight * 9949 * (lambda_**(k)) * (np.exp(-lambda_)) / math.factorial(k)



if __name__ == '__main__':
    pc_num = 1
    df_pc = pd.read_csv('D:\我的坚果云\ECS_mid\corelation.analysis\max_freq1.csv', index_col=0)
    gene_list = df_pc[df_pc['pc%d'%pc_num]>0.005].index.tolist()

    # # 混合模型分析
    # gene_list += ["0"] * (16 - len(gene_list) % 16)
    # gene_arr = np.array(gene_list).reshape(-1, 16)
    #
    # df = pd.read_csv('D:\cor_matirx\gene_filtered_result_pc%d.csv' % pc_num, index_col=0)
    # for batch in gene_arr:
    #     process = []
    #     for gene in batch:
    #         if gene != "0":
    #             data = np.array(df[gene]).reshape(-1, 1)
    #             process.append(multiprocessing.Process(target=poisson_analysis, args=(data, gene, pc_num)))
    #     for p in process:
    #         p.start()
    #     for p in process:
    #         p.join()

    # # 携带频率存在人群差异的基因
    # os.chdir('D:\我的坚果云\ECS_mid\gene_distribution_pc%d' % pc_num)
    # df_dif_gene = pd.DataFrame(columns=['average_cf', 'best_components', 'cf_list', "weights", 'highest_cf'])
    # for gene in gene_list:
    #     if gene == 'HBA1/HBA2':
    #         df = pd.read_csv('HBA1_HBA2.poisson_mixture.csv', index_col=0)
    #     else:
    #         df = pd.read_csv('%s.poisson_mixture.csv'%gene, index_col=0)
    #     best_components = df['bic'].tolist().index(df['bic'].min()) + 1
    #     df_dif_gene.loc[gene, 'average_cf'] = df.loc[1, 'parameters']
    #     df_dif_gene.loc[gene, 'best_components'] = best_components
    #     df_dif_gene.loc[gene, 'cf_list'] = df.loc[best_components, 'parameters']
    #     df_dif_gene.loc[gene, 'weights'] = df.loc[best_components, 'weights']
    #     df_dif_gene.loc[gene, 'highest_cf'] = max(df.loc[best_components, 'parameters'].split(":"))
    # df_dif_gene.to_csv('genes.different_cf.pc%d.csv'%pc_num, encoding='utf_8_sig', index=True)

    # 绘制子集中基于的携带者频数分布图
    df = pd.read_csv('D:\cor_matirx\gene_filtered_result_pc%d.csv'%pc_num, index_col=0)
    os.chdir('D:\我的坚果云\ECS_mid\gene_distribution_pc%d'%pc_num)
    df_pos = pd.read_csv('genes.different_cf.pc1.csv', index_col=0)

    for i in gene_list:
        g = sns.displot(data=df, x=i, kde=False, discrete=True)
        g.ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        par_list = df_pos.loc[i, 'cf_list'].split(":")
        weight_list = df_pos.loc[i, 'weights'].split(":")
        colors = plt.get_cmap('RdBu', len(par_list))
        for par,weight in zip(par_list, weight_list):
            y = [poisson(float(par), int(k), float(weight)) for k in np.arange(0, df[i].max()+1)]
            g.ax.plot(range(int(df[i].max()+1)), y, c=colors([par_list.index(par)]))
        plt.show()
        exit()
        if i == "HBA1/HBA2":
            plt.savefig('HBA1_HBA2.png', dpi=300)
        else:
            plt.savefig('%s.png'%i, dpi=300)
    

