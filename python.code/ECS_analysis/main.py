import copy
import os

import matplotlib.pyplot as plt
import numpy as np
import openpyxl

import data_prepare
import pandas as pd

import plot

os.environ["OMP_NUM_THREADS"] = '1'
from addresss import *


def get_risk_couples(input_df, method='any'):
    if not method in ['any', 'both', 'xlink', 'auto']:
        raise ValueError
    df = copy.deepcopy(input_df)
    df1 = df[df['member_count'] == 2]
    print("报告均已出的夫妻对数 = %d" % df1.shape[0])
    if method == 'any':
        df2 = df1[(df1['x_gene_num'] > 0) | (df1['gene_num_at_risk'] > 0)]
        return df2
    if method == 'both':
        df2 = df1[(df1['x_gene_num'] > 0) & (df1['gene_num_at_risk'] > 0)]
        return df2
    if method == 'auto':
        df2 = df1[(df1['gene_num_at_risk'] > 0)]
        return df2
    if method == 'xlink':
        df2 = df1[(df1['x_gene_num'] > 0)]
        return df2


def rm_site(input_df, gene, var_id, c_change, p_change):
    df = copy.deepcopy(input_df)
    df1 = df[df['carrier_status'] > 0]
    for i in df1[df1['var_id'].str.contains(var_id)].index:
        status_code = df.loc[i, 'carrier_status']
        auto_status, x_status, f8_inv_status, fmr1_status = 0, 0, (status_code & 2), (status_code & 1)
        gene_list = df.loc[i, 'gene'].split(":")
        gene_list.remove(gene)
        var_list = df.loc[i, 'var_id'].split(":")
        var_list.remove(var_id)
        c_change_list = df.loc[i, 'c_change'].split(":")
        c_change_list.remove(c_change)
        p_change_list = df.loc[i, 'p_change']
        p_change_list.remove(p_change)
        for t in gene_list:
            if t in Auto_list:
                auto_status = 1
            elif t in Xlink_list:
                x_status = 1
        df.loc[i, 'carrier_status'] = auto_status*2**3 + x_status*2**2 + f8_inv_status*2 + fmr1_status
        if not len(gene_list):
            df.loc[i, 'gene'] = np.nan
            df.loc[i, 'var_id'] = np.nan
            df.loc[i, 'c_change'] = np.nan
            df.loc[i, 'p_change'] = np.nan
        else:
            df.loc[i, 'gene'] = ":".join(gene_list)
            df.loc[i, 'var_id'] = ":".join(var_list)
            df.loc[i, 'c_change'] = ":".join(c_change_list)
            df.loc[i, 'p_change'] = ":".join(p_change_list)
    return df


def pos_rate_rank_by_hospital(input_df):
    df = copy.deepcopy(input_df)
    columns = ['hospital', 'reports', 'rate'] + [str(i) for i in range(1,11)]
    df2 = pd.DataFrame(columns=columns)
    hospital_list = list(set(df['hospital'].tolist()))
    gene_list = Auto_list + Xlink_list
    df_rank = pd.DataFrame(columns=['counts'], index=gene_list)
    hospital_counts = 0
    for h in hospital_list:
        for i in gene_list:
            df_rank.loc[i, 'counts'] = df[(df['hospital'] == h) & (df['gene'].str.contains(i))].shape[0]
        df_tmp = df[df['hospital'] == h]
        tmp_list = [df_tmp[df_tmp['carrier_status'] > 0].shape[0], df_tmp.shape[0]]
        tmp_list[1] = tmp_list[0] / tmp_list[1]
        df2.loc[hospital_counts] = [h] + tmp_list + df_rank.sort_values('counts', ascending=False).index.tolist()[0:10]
        hospital_counts += 1
    return df2


if __name__ == '__main__':
    # 分别读取男女数据，合并为个体为单位的表格
    # df_id_table = pd.read_excel("primary.data/国家课题身份证前六位-截至20230403.xlsx")
    # df11 = pd.read_excel("primary.data/项目编号_Y0013导出1.xlsx")
    # df12 = pd.read_excel("primary.data/项目编号_Y0013导出2.xlsx")
    # df21 = pd.read_excel("primary.data/项目编号_Y0014导出1.xlsx")
    # df22 = pd.read_excel("primary.data/项目编号_Y0014导出2.xlsx")
    # df11 = data_prepare.convert_in_samples(df11, sex_label=0, id_table=df_id_table)
    # df12 = data_prepare.convert_in_samples(df12, sex_label=0, id_table=df_id_table)
    # df21 = data_prepare.convert_in_samples(df21, sex_label=1, id_table=df_id_table)
    # df22 = data_prepare.convert_in_samples(df22, sex_label=1, id_table=df_id_table)
    # df = pd.concat([df11, df12, df21, df22])
    # df.drop_duplicates(inplace=True)
    # df = df[~df['name'].str.contains("测试")]
    # df.to_csv("sample.csv", index=False, encoding='utf_8_sig')

    # df_area = pd.read_csv("area.csv", index_col='area')
    # df_area = data_prepare.transform_merge_area(df_area, Area_counterparts)
    # df_area.drop('unknown', inplace=True)
    # plot.plot_gene(df_area)
    # plot.plot_area_individual(df_area)
    # plot.plot_area_pca(df_area)
    # plot.plot_kmeans_pca(df_area, 3)
    # plot.kmeans_evaluations(df_area)
    # plot.plot_area_area_heatmap(df_area)
    # plot.plot_area2_fst_heatmap(df_area)
    # plot.plot_ckmeans(df_area, 7)
    # plot.plot_ckmeans_evaluation(df_area)

    df = pd.read_csv('sample.csv')
    # df = df[~df['c_change'].isna()]
    # gene_list = []
    # var_list = []
    # c_change_list = []
    # p_change_list = []
    # for i in df.index:
    #     var_list += df.loc[i, 'var_id'].split(":")
    #     gene_list += df.loc[i, 'gene'].split(":")
    #     c_change_list += df.loc[i, 'c_change'].split(":")
    #     p_change_list += df.loc[i, 'p_change'].split(":")
    # gene_list_rmdup = list(set(gene_list))
    # gene_list_counts = [gene_list.count(i) for i in gene_list_rmdup]
    # var_list_rmdup = list(set(var_list))
    # var_list_counts = [var_list.count(i) for i in var_list_rmdup]
    # df_count = pd.DataFrame(columns=['var_id', 'gene', 'counts', 'all_var_in_gene', 'ratio', 'c_change', 'p_change'])
    # for i in range(len(var_list_rmdup)):
    #     index_ = var_list.index(var_list_rmdup[i])
    #     df_count.loc[i, 'var_id'] = var_list_rmdup[i]
    #     df_count.loc[i, 'gene'] = gene_list[index_]
    #     df_count.loc[i, 'c_change'] = c_change_list[index_]
    #     df_count.loc[i, 'p_change'] = p_change_list[index_]
    #     df_count.loc[i, 'counts'] = var_list_counts[i]
    #     df_count.loc[i, 'all_var_in_gene'] = gene_list_counts[gene_list_rmdup.index(gene_list[index_])]
    #
    # df_count.to_csv('tmp.csv', encoding='utf_8_sig', index=False)
    # print(df_count)
    print(df)

