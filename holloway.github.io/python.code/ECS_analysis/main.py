import os
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
    # df = pd.read_excel('tmp.xlsx', index_col=0)
    #
    # fig = plt.figure(figsize=(8,4))
    # sns.set_theme(style='ticks')
    # ax = fig.add_axes((0.15, 0.3, 0.3, 0.5))
    # sns.heatmap(data=df, cmap='coolwarm', linewidths=1, )
    # ax.set_xlabel('Monthly Household Income(CHY)', fontsize=10)
    # ax.set_ylabel('Severity of Phenotype', fontsize=10)
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=10)
    # ax.set_yticklabels(ax.get_yticklabels(), fontsize=10)
    # ax.collections[0].colorbar.ax.tick_params(labelsize=10)
    # ax.collections[0].colorbar.set_label('Support Rate for Inclusion', fontsize=10)
    # plt.show()

    # df = pd.read_excel('demographics.xlsx', sheet_name=None)
    x = np.array([1, 2, 3, 4, 5])  # x坐标
    y = np.array([2.3, 3.1, 4.5, 3.8, 5.2])  # y坐标
    y_err = np.array([0.1, 0.2, 0.15, 0.25, 0.18])  # y坐标的误差值
    fig, ax = plt.subplots()
    ax.errorbar(x, y, yerr=y_err, fmt='o-', capsize=5)
    ax.set_title('Line Plot with Error Bars')
    ax.set_xlabel('X Axis Label')
    ax.set_ylabel('Y Axis Label')
    plt.show()
