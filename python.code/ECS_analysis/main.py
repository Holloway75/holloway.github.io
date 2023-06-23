import copy
import os
import numpy as np
from pyckmeans import CKmeans
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
import data_prepare
import plot
from addresss import *
import pandas as pd
from PoissonMixture import PoissonMixture
import multiprocessing


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
            raise ValueError
    return df


if __name__ == '__main__':
    os.chdir('D:\我的坚果云\ECS_1.6w_samples\gene_distribution_pc1')
    df = pd.read_csv('D:\我的坚果云\ECS_1.6w_samples\gene_filtered_result_pc1.csv', index_col=0)
    # gene_list = []
    #
    # for gene in df.columns:
    #     if df[gene].max() > 5:
    #         gene_list.append(gene)
    # print(len(gene_list))

    # # 混合模型分析
    # gene_list += ["0"] * (16 - len(gene_list) % 16)
    # gene_arr = np.array(gene_list).reshape(-1, 16)
    #
    # for batch in gene_arr:
    #     process = []
    #     for gene in batch:
    #         if gene != "0":
    #             data = np.array(df[gene]).reshape(-1, 1)
    #             process.append(multiprocessing.Process(target=poisson_mixture_analysis, args=(data, gene)))
    #     for p in process:
    #         p.start()
    #     for p in process:
    #         p.join()

    # # 携带频率存在人群差异的基因
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
    # df_dif_gene.to_csv('genes.different_cf.pc1.csv', encoding='utf_8_sig', index=True)

    # 绘制子集中携带者频数分布图
    # df_pos = pd.read_csv('genes.different_cf.pc1.csv', index_col=0)
    # for i in gene_list:
    #     g = sns.displot(data=df, x=i, kde=False, discrete=True)
    #     g.ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    #
    #     par_list = df_pos.loc[i, 'cf_list'].split(":")
    #     weight_list = df_pos.loc[i, 'weights'].split(":")
    #     colors = plt.get_cmap('RdBu', len(par_list))
    #     for par,weight in zip(par_list, weight_list):
    #         y = [poisson(float(par), int(k), float(weight)) for k in np.arange(0, df[i].max()+1)]
    #         g.ax.plot(np.arange(0, df[i].max()+1), y, c=plt.cm.Accent(par_list.index(par)), label='lambda=%.2f weight=%.2f'%(float(par), float(weight)))
    #     plt.legend()
    #     if i == "HBA1/HBA2":
    #         plt.savefig('HBA1_HBA2.png', dpi=300)
    #     else:
    #         plt.savefig('%s.png'%i, dpi=300)

    # gene_list = ['USH2A', 'GAA', 'TYR']
    # fig, ax = plt.subplots(1, 3, constrained_layout=True, figsize=(15, 5))
    #
    # for i in gene_list:
    #     axesSub = sns.histplot(data=df, x=i, ax=ax[gene_list.index(i)], discrete=True)
    #     par_list = df_pos.loc[i, 'cf_list'].split(":")
    #     weight_list = df_pos.loc[i, 'weights'].split(":")
    #     colors = plt.get_cmap('RdBu', len(par_list))
    #     lines = []
    #     for par,weight in zip(par_list, weight_list):
    #         y = [poisson(float(par), int(k), float(weight)) for k in np.arange(0, df[i].max()+1)]
    #         line, = ax[gene_list.index(i)].plot(np.arange(0, df[i].max()+1), y, c=plt.cm.Accent(par_list.index(par)), label='lambda=%.2f weight=%.2f'%(float(par), float(weight)))
    #         lines.append(line)
    # for a in ax:
    #     a.legend()
    # plt.show()

    path = r'D:\ECS报告\20230621'
    df2 = get_sample_name(path)
    print(df2)
