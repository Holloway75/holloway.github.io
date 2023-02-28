import pandas as pd

from addresss import *
from plot import *


def get_areas_from_id(name, fid, data):
    df1 = data[data.fid == fid]
    if not df1.shape[0]:
        return "unknown"
    df2 = df1[df1.name == name]
    if not df2.shape[0]:
        return "unknown"
    df2.index = [0]
    if pd.isna(df2.at[0, 'id']):
        return "unknown"
    else:
        b = str(df2.at[0, 'id'])
        return id_province[b[0:2]]


def convert_in_samples(df, data, alist, xlist, sex_label):
    """
    将原始df数据转录到新df2, df每行为一个个体
    :param df: 原始数据
    :param df2: 新数据,每个个体为一行
    :param data: 包含name, fid与身份证前六位的表
    :param alist: 常染色体基因列表
    :param xlist: x染色体基因列表
    """
    columns_sample = ['name', 'sex', 'fid', 'area', 'hospital', 'carrier_status', 'gene', 'var_id']
    df2 = pd.DataFrame(columns=columns_sample)
    sample_counts = 0
    lines = 0
    all_lines = df.shape[0]
    while lines < all_lines:
        t = 1
        if not lines == all_lines-1:
            while pd.isna(df.loc[lines+t, '姓名']):
                t += 1
                if lines + t == all_lines:
                    break
        df2.loc[sample_counts, 'name'] = df.at[lines, '姓名']
        df2.loc[sample_counts, 'sex'] = sex_label
        df2.loc[sample_counts, 'fid'] = df.at[lines, '家系编号']
        df2.loc[sample_counts, 'area'] = get_areas_from_id(df.loc[lines, '姓名'], df.at[lines, '家系编号'], data)
        df2.loc[sample_counts, 'hospital'] = df.at[lines, '送检医院']

        auto_status, x_status, f8_inv_status, fmr1_status = 0,0,0,0
        sample_gene_list, sample_var_list = [], []

        if df.loc[lines, '内含子1'] == '杂合变异' or df.loc[lines, '内含子22'] == '杂合变异':
            f8_inv_status = 1
            sample_gene_list.append("F8")
            sample_var_list.append("F8_inv")

        if not pd.isna(df.loc[lines, 'CGG重复数目']):
            cgg_num = df.loc[lines, 'CGG重复数目'].split('|')
            if int(cgg_num[0]) > 54 or int(cgg_num[1]) > 54:
                fmr1_status = 1
                sample_gene_list.append("FMR1")
                sample_var_list.append(df.loc[lines, 'CGG重复数目'])

        if not pd.isna(df.loc[lines, '基因']):
            for i in range(t):
                sample_gene_list.append(df.loc[lines+i, '基因'])
                sample_var_list.append(df.loc[lines+i, '变异ID'])

        if sum([(x in alist) for x in sample_gene_list]):
            auto_status = 1
        if sum([(x in xlist) for x in sample_gene_list]):
            x_status = 1
        carrier_status = auto_status*2**3 + x_status*2**2 + f8_inv_status*2 + fmr1_status
        df2.loc[sample_counts, 'carrier_status'] = carrier_status
        df2.loc[sample_counts, 'gene'] = ":".join(sample_gene_list)
        df2.loc[sample_counts, 'var_id'] = ":".join(sample_var_list)
        lines += t
        sample_counts += 1
    return df2





if __name__ == '__main__':
    # df_sample = pd.read_csv('sample.combined.csv')
    # df = pd.DataFrame(columns=[str(i) for i in range(100)])
    # for i in range(150, 210, 10):
    #     df.loc[i] = [data_prepare.random_n_distance(i, df_sample) for t in range(100)]
    # df.to_csv('n_test1.csv', index=True)

    df = pd.read_csv('n.test.10_210_10.csv', index_col=0)
    g = sns.catplot(
        data=df.T, kind="strip"
    )
    plt.show()






