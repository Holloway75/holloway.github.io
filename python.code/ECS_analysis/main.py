import copy
import os
import numpy as np
from addresss import *
import pandas as pd
import addresss
import data_prepare
import plot
import seaborn as sns
from matplotlib import pyplot as plt


def to_compare(input_df_, pc_num):
    input_df = copy.deepcopy(input_df_)
    sample_index_list = input_df.index.tolist()
    pc_list = ['PC%d' % i for i in range(1, pc_num+1)]
    df_cor_matrx = pd.DataFrame(columns=sample_index_list)
    df_gene_filtered_result = pd.DataFrame(index=Auto_list, columns=sample_index_list, dtype='int8')
    for arr in sample_index_list:
        arr_loc = np.array(input_df.loc[arr, pc_list])
        for col in sample_index_list:
            col_loc = np.array(input_df.loc[col, pc_list])
            df_cor_matrx.loc[arr, col] = np.linalg.norm(arr_loc-col_loc)
        nearest_1000_list = df_cor_matrx.loc[arr].sort_values().index.tolist()[0:1000]
        input_df.loc[arr, 'nearest_1000_list'] = ":".join(nearest_1000_list)
        df_tmp = input_df[[(i in nearest_1000_list) for i in input_df.index]]
        gene_list = df_tmp[~df_tmp['gene'].isna()]['gene'].tolist()
        df_gene_filtered_result[arr] = [gene_list.count(i) for i in Auto_list]
    df_cor_matrx.to_csv('cor_matrix_pc%d.csv'%pc_num, index=True, encoding='utf_8_sig')
    input_df.to_csv('sample.nearest_1000_by-pc%d.csv' % pc_num, index=True, encoding='utf_8_sig')
    df_gene_filtered_result.to_csv('gene_filtered_result_pc%d.csv', index=True, encoding='utf_8_sig')


if __name__ == '__main__':
    os.chdir('D:\我的坚果云\ECS_mid')
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

    df_sample = pd.read_csv('sample.coord.csv', index_col='indivID')
    df_sample.drop(index=df_sample[df_sample['main_area']=='unknown'].index, inplace=True)
    to_compare(df_sample, 20)
    # print(df_sample.loc[0, ['PC1','PC2']])

