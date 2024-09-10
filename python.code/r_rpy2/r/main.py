import numpy as np
import os
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr
import statsmodels.api as sm
stats = importr('stats')

if __name__ == '__main__':
    os.chdir('E:\我的坚果云\胎盘炎症-极早产')
    df_x1 = pd.read_excel('primary_data.xlsx', index_col=0, sheet_name='病理特征')
    df_x2 = pd.read_excel('primary_data.xlsx', index_col=0, sheet_name='一般特征')
    df_y = pd.read_excel('primary_data.xlsx', index_col=0, sheet_name='结局')

    input_data = [df_x1, df_x2, df_y]

    def merge_twins(df_list):
        _df_x1 = pd.read_excel('primary_data.xlsx', index_col=0, sheet_name='双胎病理')
        _df_x2 = pd.read_excel('primary_data.xlsx', index_col=0, sheet_name='双胎一般特征')
        _df_y = pd.read_excel('primary_data.xlsx', index_col=0, sheet_name='双胎结局')

        _df_x1 = pd.concat([df_list[0], _df_x1], axis=0)
        _df_x2 = pd.concat([df_list[1], _df_x2], axis=0)
        _df_y = pd.concat([df_list[2], _df_y], axis=0)

        return [_df_x1, _df_x2, _df_y]

    def data_preprocess(df_list):
        for _df in df_list:
            _df.fillna(value=0, inplace=True)
            _df.replace({'无': 0, '有': 1}, inplace=True)

    input_data = merge_twins(input_data)
    data_preprocess(input_data)

    df_x = pd.merge(input_data[0], input_data[1], left_index=True, right_index=True)
    y_ = input_data[2]['支气管肺发育不良']

    col_for_std = ['产妇分娩时年龄', '分娩孕周']

    df_stats = pd.DataFrame(index=df_x.columns, columns=['statistic', 'p', 'method'])
    for col in df_x.columns:
        x_ = df_x[col]
        if col in col_for_std:
            df_stats.loc[col, 'method'] = 'Logit'
            x_ = sm.add_constant(x_)
            res = sm.Logit(y_, x_).fit()
            df_stats.loc[col, 'statistic'] = res.tvalues[col]
            df_stats.loc[col, 'p'] = res.pvalues[col]
        else:
            cross = pd.crosstab(y_, x_)
            res = stats.fisher_test(cross)
            df_stats.loc[col, 'method'] = 'fisher_test'
            df_stats.loc[col, 'p'] = res[0][0]

    df_stats.to_excel('stats1.xlsx', index=True)