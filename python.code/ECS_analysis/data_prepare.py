from addresss import *
import fst_calculation as fst


def data2plot_gene(input_df, cut_line=1/200, area=None):
    if not area is None:
        pre_df = input_df.loc[area]
    else:
        pre_df = copy.deepcopy(input_df)
    pre_df = fst.filter_by_cf(pre_df, cut_line)

    # 增加一行'total' 并统计总人数
    pre_df.loc['total'] = [sum(pre_df[t]) for t in pre_df.columns]
    male_counts = pre_df.loc['total', 'individuals_male']
    total_counts = pre_df.loc['total', 'individuals_total']
    female_counts = total_counts - male_counts
    glist = pre_df.columns.tolist()[5:]

    # 常隐和x连锁分别计算携带频率
    for gene in glist:
        if gene in auto_list:
            pre_df.loc['total', gene] = pre_df.loc['total', gene] / total_counts
        elif gene in xlink_list:
            pre_df.loc['total', gene] = pre_df.loc['total', gene] / female_counts
        else:
            raise ValueError

    # 去除人数相关的列，只保留基因，并以携带频率排序
    pre_df = pre_df[glist]
    pre_df.sort_values(by='total', inplace=True, axis=1)
    return pre_df

