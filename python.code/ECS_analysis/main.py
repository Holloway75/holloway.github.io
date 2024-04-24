import os

import pandas

import data_prepare
import plot
os.environ["OMP_NUM_THREADS"] = '1'
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from addresss import *
import pandas as pd
from sklearn.preprocessing import OneHotEncoder
import statsmodels.api as sm

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
        df.to_csv('E:\我的坚果云\ECS_1.6w_samples_corrected\gene_distribution_pc2_auto\HBA1_HBA2.poisson_mixture.csv', encoding='utf_8_sig')
    else:
        df.to_csv('E:\我的坚果云\ECS_1.6w_samples_corrected\gene_distribution_pc2_auto\%s.poisson_mixture.csv' % gene, encoding='utf_8_sig')


def poisson(lambda_, k, weight):
    lambda_, k, weight = float(lambda_), float(k), float(weight)
    log_prob = -lambda_ + k * np.log(lambda_) - sum([np.log(i) for i in np.arange(1, k+1)])
    return weight * 15089 * np.exp(log_prob)


def fac(n):
    return np.log(np.arange(1, n+1)).sum()


def binominal(p, k, weight):
    p, k, weight = float(p), float(k), float(weight)
    log_prob = fac(1000) - fac(k) - fac(1000-k) + k * np.log(p) + (1000-k) * np.log(1-p)
    return weight * 15089 * np.exp(log_prob)


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


def get_risk(varid, table:pd.DataFrame):
    gene = table.loc[varid, 'gene']
    var_count = table.loc[varid, 'counts']
    df_tmp = table[table.gene == gene]
    if var_count == df_tmp['counts'].max:
        car_count_pre = 0
    else:
        car_count_pre = df_tmp[df_tmp.counts > var_count]['counts'].sum()
    car_count_after = car_count_pre + var_count

    if gene in Auto_list:
        return (car_count_after/33086)**2 - (car_count_pre/33086)**2
    elif gene in Xlink_list:
        return car_count_after/16640 - car_count_pre/16640
    else:
        raise ValueError


if __name__ == '__main__':
    os.chdir('E:\我的坚果云\ecs问卷')
    df = pd.read_excel('tmp.xlsx', index_col=0)

    fig = plt.figure(figsize=(8,4))
    sns.set_theme(style='ticks')
    # fig1
    ax = fig.add_axes((0.15, 0.3, 0.3, 0.55))
    sns.heatmap(data=df, cmap='coolwarm', linewidths=1, )
    ax.set_xlabel('Monthly Household Income(CHY)', fontsize=10)
    ax.set_ylabel('Severity of Phenotype', fontsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)
    ax.collections[0].colorbar.ax.tick_params(labelsize=10)
    ax.collections[0].colorbar.set_label('Support Rate for Inclusion', fontsize=10)
    ax.set_title('a', x=0.5, y=-0.46, fontsize=10)

    # fig2
    df = pd.read_excel('demographics.xlsx', sheet_name=None)
    x = np.arange(1, 6)
    y_1 = []
    y_1_err = [[], []]
    y_2 = []
    y_2_err = [[], []]
    y_3 = []
    y_3_err = [[], []]
    y_4 = []
    y_4_err = [[], []]
    tar = [y_1, y_1_err, y_2, y_2_err, y_3, y_3_err, y_4, y_4_err]

    def add_value(df:pandas.DataFrame, l1, l1_err, l2, l2_err, l3, l3_err, l4, l4_err):
        l1.append(df.loc['5000元以下', 'OR值'])
        l1_err[0].append(df.loc['5000元以下', 'lower'])
        l1_err[1].append(df.loc['5000元以下', 'upper'])
        l2.append(df.loc['5000元-1万元', 'OR值'])
        l2_err[0].append(df.loc['5000元-1万元', 'lower'])
        l2_err[1].append(df.loc['5000元-1万元', 'upper'])
        l3.append(df.loc['2万-3万', 'OR值'])
        l3_err[0].append(df.loc['2万-3万', 'lower'])
        l3_err[1].append(df.loc['2万-3万', 'upper'])
        l4.append(df.loc['3万以上', 'OR值'])
        l4_err[0].append(df.loc['3万以上', 'lower'])
        l4_err[1].append(df.loc['3万以上', 'upper'])


    ax = fig.add_axes((0.58, 0.22, 0.37, 0.7))
    pal = sns.color_palette(n_colors=4)
    x_ = np.array([1, 2, 3, 4, 5])  # x坐标
    for i in np.arange(1, 6):
        df_tmp = df[list(df.keys())[i]].set_index('社会人口因素')
        add_value(df_tmp, y_1, y_1_err, y_2, y_2_err, y_3, y_3_err, y_4, y_4_err)

    labels = ['<5k', '5k-10k', '20k-30k', '>30k', '10k-20k\n(contrast)']
    for i in range(4):
        tar[2*i+1][0] = [a-b for a, b in zip(tar[2*i], tar[2*i+1][0])]
        tar[2 * i + 1][1] = [a-b for a, b in zip(tar[2 * i + 1][1], tar[2*i])]
        ax.errorbar(x_+(i-2)/7, tar[2*i], yerr=tar[2*i+1], fmt='o-', capsize=5, c=pal[i])
        ax.plot(x_+(i-2)/7, tar[2*i], c=pal[i], label=labels[i])

    ax.hlines(y=1, xmin=0.6, xmax=5.5, label=labels[4], colors='black')
    ax.set_xlim(0.5, 5.5)
    ax.set_ylim(0.5, 2)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xlabel('Severity of Phenotype', fontsize=10)
    ax.set_xticks(np.arange(1, 6))
    ax.set_xticklabels([i.capitalize() for i in list(df.keys())[1:]], fontsize=10)
    ax.set_ylabel('OR', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)
    legend = ax.legend(title='Monthly Household Income(CHY)', fontsize=10, labelspacing=0.5, ncol=2,
                       columnspacing=0.5, handletextpad=0.4)
    legend.get_title().set_fontsize(fontsize=10)
    ax.set_title('b', x=0.5, y=-0.25, fontsize=10)
    plt.show()
