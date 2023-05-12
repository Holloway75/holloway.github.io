import os
import numpy as np
from addresss import *
import pandas as pd
import addresss
import data_prepare
import plot
import seaborn as sns
from matplotlib import pyplot as plt


def to_compare(input_df, pc_num):
    sample_index_list = input_df.index.tolist()
    pc_list = ['PC%d' % i for i in range(1, pc_num+1)]
    df_cor_matrx = pd.DataFrame(columns=sample_index_list)
    for arr in sample_index_list:
        print(arr)
        arr_loc = np.array(input_df.loc[arr, pc_list])
        for col in sample_index_list:
            col_loc = np.array(input_df.loc[col, pc_list])
            df_cor_matrx.loc[arr, col] = np.linalg.norm(arr_loc-col_loc)
        nearest_1000_list = df_cor_matrx.loc[arr].sort_values().index.tolist()[0:1000]
        input_df.loc[arr, 'nearest_1000_list'] = ":".join(nearest_1000_list)
        print(input_df)
        exit()



if __name__ == '__main__':
    os.chdir('D:\我的坚果云\ECS_mid')
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

    df_sample = pd.read_csv('sample.coord.csv', index_col='indivID')
    df_sample.drop(index=df_sample[df_sample['main_area']=='unknown'].index, inplace=True)
    to_compare(df_sample, 20)
    # print(df_sample.loc[0, ['PC1','PC2']])

