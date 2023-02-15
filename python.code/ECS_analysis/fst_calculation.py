<<<<<<< HEAD
=======
import copy

>>>>>>> f85d5d3f0f61e82a369e5026d6023a7710853ed9
from addresss import *


def get_fst_from_area2_gene(input_df, area_arr, area_col, gene):
    if gene in auto_list:
        samples_arr = input_df.loc[area_arr, 'individuals_total']
        samples_col = input_df.loc[area_col, 'individuals_total']
    elif gene in xlink_list:
        samples_arr = input_df.loc[area_arr, 'individuals_total'] - input_df.loc[area_arr, 'individuals_male']
        samples_col = input_df.loc[area_col, 'individuals_total'] - input_df.loc[area_col, 'individuals_male']
    else:
        raise ValueError

    freq_arr = input_df.loc[area_arr, gene] / (2*samples_arr)
    freq_col = input_df.loc[area_col, gene] / (2*samples_col)
    # hs = (2*freq_arr*(1-freq_arr)*samples_arr + 2*freq_col*(1-freq_col)*samples_col) / (samples_arr + samples_col)
    # freq_total = (input_df.loc[area_arr, gene] + input_df.loc[area_col, gene]) / (2*(samples_arr + samples_col))

    hs = freq_arr*(1-freq_arr) + freq_col*(1-freq_col)
    freq_total = (freq_arr + freq_col)/2
    ht = 2*freq_total*(1-freq_total)
    if not ht:
        return 0
    fst = 1-hs/ht
    if fst < 10**-12:
        return 0
    else:
        return fst


def get_average_fst_from_area2(input_df, area_arr, area_col):
    a = 0
    glist = input_df.columns.tolist()[5:]
    for gene in glist:
        a += get_fst_from_area2_gene(input_df, area_arr, area_col, gene)
    return a/len(glist)


def filter_by_cf(input_df, cut_line=1/200):
    pre_df = copy.deepcopy(input_df)
    rm_list = []
    for gene in auto_list:
        cf = sum(input_df[gene].tolist()) / sum(input_df['individuals_total'])
        if cf < cut_line:
            rm_list.append(gene)
    for gene in xlink_list:
        cf = sum(input_df[gene].tolist()) / (sum(input_df['individuals_total']) - sum(input_df['individuals_male']))
        if cf < cut_line:
            rm_list.append(gene)
    pre_df.drop(rm_list, inplace=True, axis=1)
    return pre_df


def data_prepare_for_heatmap(pre_df):
    input_df = copy.deepcopy(pre_df)
    input_df = filter_by_cf(input_df, cut_line=1/200)
    area_heatmap_list = input_df.index.tolist()
    pre_df2 = pd.DataFrame()
    for arr in area_heatmap_list:
        for col in area_heatmap_list:
            pre_df2.loc[arr, col] = get_average_fst_from_area2(input_df, arr, col)
    pre_df2.index = pre_df2.index.astype('category').set_categories(area_sort_list, ordered=True)
    pre_df2.sort_index(inplace=True)
    return pre_df2[area_sort_list]

