import numpy as np
from addresss import *
import copy
import pandas as pd


def get_fst_from_area2_gene(input_df, area_row, area_col, gene):
    if gene in Auto_list:
        samples_row = input_df.loc[area_row, 'individuals_total']
        samples_col = input_df.loc[area_col, 'individuals_total']
    elif gene in Xlink_list:
        samples_row = input_df.loc[area_row, 'individuals_total'] - input_df.loc[area_row, 'individuals_male']
        samples_col = input_df.loc[area_col, 'individuals_total'] - input_df.loc[area_col, 'individuals_male']
    else:
        raise ValueError

    cf_row = input_df.loc[area_row, gene] / (samples_row)
    cf_col = input_df.loc[area_col, gene] / (samples_col)
    af_row = cf_row/2
    af_col = cf_col/2
    af_total = 0.5*(af_col + af_row)

    hs = af_row * (1-af_row) + af_col * (1-af_col)
    ht = 2 * af_total * (1-af_total)

    if not ht:
        return np.nan
    fst = 1-hs/ht
    if fst < 10**-12:
        return 0
    else:
        return fst


def get_average_fst_from_area2(input_df, area_arr, area_col):
    a = []
    glist = input_df.columns.tolist()[5:]
    for gene in glist:
        a = np.append(a, get_fst_from_area2_gene(input_df, area_arr, area_col, gene))
        a = np.array(a)
        a = a[~np.isnan(a)]
    return a.mean()


def filter_by_cf(input_df, cut_line=1/200):
    pre_df = copy.deepcopy(input_df)
    rm_list = []
    glist = pre_df.columns.tolist()[5:]

    # 常隐和x连锁分别计算携带频率
    for gene in glist:
        if gene in Auto_list:
            cf = sum(input_df[gene].tolist()) / sum(input_df['individuals_total'])
            if cf < cut_line:
                rm_list.append(gene)
        elif gene in Xlink_list:
            cf = sum(input_df[gene].tolist()) / (sum(input_df['individuals_total']) - sum(input_df['individuals_male']))
            if cf < cut_line:
                rm_list.append(gene)
        else:
            raise ValueError

    pre_df.drop(rm_list, inplace=True, axis=1)
    return pre_df


def data_prepare_for_heatmap(pre_df, cut_line=1/200):
    input_df = copy.deepcopy(pre_df)
    input_df = filter_by_cf(input_df, cut_line)
    area_heatmap_list = input_df.index.tolist()
    pre_df2 = pd.DataFrame()
    for arr in area_heatmap_list:
        for col in area_heatmap_list:
            pre_df2.loc[arr, col] = get_average_fst_from_area2(input_df, arr, col)
    pre_df2.sort_values(by='陕甘宁', axis=1, inplace=True)
    pre_df2.sort_values(by='陕甘宁', axis=0, inplace=True)
    return pre_df2

