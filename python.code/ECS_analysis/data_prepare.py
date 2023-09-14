from addresss import *
import fst_calculation as fst
import pandas as pd
import copy
import numpy as np
import re


__all__ = [
    'convert_in_samples'
]


def get_keys(dic, val):
    for k, v in dic.items():
        if val in v:
            return k
    return val


def convert_in_samples(input_df, sex_label):
    """
    项目导出文件转换格式，文件内样本性别相同，
    :param input_df: 项目导出文件，0013为女性，0014为男性，每个样本占用行数为检出变异数，至少1行，单个样本可能占用多行
    :param sex_label: 性别标签，男性为'1'，女性为'0'
    :param id_table: 包含姓名、家系id、身份证前六位的表
    :return:转化后格式，每个样本占1行
    """
    columns_sample = ['name', 'sex', 'fid', 'area', 'hospital', 'carrier_status', 'gene', 'var_id', 'c_change',
                      'p_change']
    df2 = pd.DataFrame(columns=columns_sample)
    sample_counts = 0
    all_lines = input_df.shape[0]
    index = 0
    while index < all_lines:
        t = 1  # t为当前样本所占用的行数
        if not index == all_lines - 1:
            while pd.isna(input_df.loc[index + t, '姓名']):
                t += 1
                if index + t == all_lines:
                    break
        df2.loc[sample_counts, 'name'] = input_df.loc[index, '姓名']
        df2.loc[sample_counts, 'sex'] = sex_label
        df2.loc[sample_counts, 'fid'] = input_df.loc[index, '家系编号']
        df2.loc[sample_counts, 'hospital'] = input_df.loc[index, '送检医院']

        auto_status, x_status, f8_inv_status, fmr1_status = 0, 0, 0, 0
        sample_gene_list, sample_var_list, c_change_list, p_change_list = [], [], [], []

        if not sex_label:
            if input_df.loc[index, '内含子1'] == '杂合变异' or input_df.loc[index, '内含子22'] == '杂合变异':
                f8_inv_status = 1
                sample_gene_list.append("F8")
                sample_var_list.append("F8_inv")
                c_change_list.append("F8_inv")
                p_change_list.append("F8_inv")

            if not pd.isna(input_df.loc[index, 'CGG重复数目']):
                cgg_max = max(list(map(int, list(filter(None, re.split('[^0-9]', input_df.loc[index, 'CGG重复数目']))))))
                if cgg_max > 54:
                    fmr1_status = 1
                    sample_gene_list.append("FMR1")
                    sample_var_list.append(input_df.loc[index, 'CGG重复数目'])
                    c_change_list.append("FMR1_CGG重复数目%d" % cgg_max)
                    p_change_list.append("FMR1_CGG重复数目%d" % cgg_max)

        if not pd.isna(input_df.loc[index, '基因']):
            for i in range(t):
                sample_gene_list.append(input_df.loc[index+i, '基因'])
                sample_var_list.append(input_df.loc[index+i, '变异ID'])
                c_change_list.append(input_df.loc[index+i, '核苷酸改变'])
                p_change_list.append(input_df.loc[index+i, '氨基酸改变'])

        if sum([(x in Auto_list) for x in sample_gene_list]):
            auto_status = 1
        if sum([(x in Xlink_list) for x in sample_gene_list]):
            x_status = 1
        carrier_status = auto_status*2**3 + x_status*2**2 + f8_inv_status*2 + fmr1_status
        df2.loc[sample_counts, 'carrier_status'] = carrier_status
        df2.loc[sample_counts, 'gene'] = ":".join(sample_gene_list)
        df2.loc[sample_counts, 'var_id'] = ":".join(sample_var_list)
        df2.loc[sample_counts, 'c_change'] = ":".join(c_change_list)
        df2.loc[sample_counts, 'p_change'] = ":".join(p_change_list)
        index += t
        sample_counts += 1
    return df2


def convert_in_couples(input_df):
    columns = ['fid', 'member_count', 'gene_at_risk', 'gene_num_at_risk', 'female_name', 'female_carrier_status',
               'female_gene', 'female_variant', 'male_name', 'male_carrier_status', 'male_gene', 'male_variant',
               'x_gene', 'x_gene_num', 'hospital']
    fid_list =list(set(input_df['fid'].tolist()))
    df_couple = pd.DataFrame(columns=columns)

    fam_count = 0
    for i in fid_list:
        df_tmp = input_df[input_df['fid'] == i]

        df_couple.loc[fam_count, 'fid'] = i
        df_couple.loc[fam_count, 'member_count'] = df_tmp.shape[0]
        if df_tmp.shape[0] > 2:
            print(df_tmp)
            raise ValueError

        for t in df_tmp.index:
            if df_tmp.loc[t, 'sex']:
                df_couple.loc[fam_count, 'male_name'] = df_tmp.loc[t, 'name']
                df_couple.loc[fam_count, 'male_carrier_status'] = df_tmp.loc[t, 'carrier_status']
                df_couple.loc[fam_count, 'male_gene'] = df_tmp.loc[t, 'gene']
                df_couple.loc[fam_count, 'male_variant'] = df_tmp.loc[t, 'var_id']
            else:
                df_couple.loc[fam_count, 'female_name'] = df_tmp.loc[t, 'name']
                df_couple.loc[fam_count, 'female_carrier_status'] = df_tmp.loc[t, 'carrier_status']
                df_couple.loc[fam_count, 'female_gene'] = df_tmp.loc[t, 'gene']
                df_couple.loc[fam_count, 'female_variant'] = df_tmp.loc[t, 'var_id']
        if df_couple.loc[fam_count, 'member_count'] == 2:
            if not pd.isna(df_couple.loc[fam_count, 'female_gene']):
                female_gene_list = df_couple.loc[fam_count, 'female_gene'].split(':')
                x_gene = list(set(female_gene_list) & set(Xlink_list))
                df_couple.loc[fam_count, 'x_gene'] = ":".join(x_gene)
                df_couple.loc[fam_count, 'x_gene_num'] = len(x_gene)
                if not pd.isna(df_couple.loc[fam_count, 'male_gene']):
                    male_gene_list = df_couple.loc[fam_count, 'male_gene'].split(':')
                    gene_risk_list = list(set(female_gene_list) & set(male_gene_list))
                    df_couple.loc[fam_count, 'gene_at_risk'] = ":".join(gene_risk_list)
                    df_couple.loc[fam_count, 'gene_num_at_risk'] = len(gene_risk_list)
        df_couple.loc[fam_count, 'hospital'] = df_tmp['hospital'].tolist()[0]
        fam_count += 1
    return df_couple


def convert_in_areas(df):
    constant_column = ['carriers_auto', 'carriers_x', 'carriers_total',
                       'individuals_male', 'individuals_total']
    title_gene_list = Auto_list + Xlink_list
    column = constant_column + title_gene_list
    areas = list(set(df.area.tolist()))
    df2 = pd.DataFrame(index=areas, columns=column)
    df2.index.name = 'area'
    for area in df2.index:
        df_tmp = df[df.area == area]

        # 计数该地区中 常染色体病携带者 x携带者 总携带者人数 男性人数 总人数
        status_list = df_tmp.carrier_status.tolist()
        df2.loc[area, 'carriers_auto'] = len([sam for sam in status_list if (sam & 8)])
        df2.loc[area, 'carriers_x'] = len([sam for sam in status_list if (sam & 4)])
        df2.loc[area, 'carriers_total'] = len([sam for sam in status_list if sam])
        df2.loc[area, 'individuals_male'] = df_tmp.sex.sum()
        df2.loc[area, 'individuals_total'] = len(status_list)

        # 计数该地区中每个基因的出现次数
        df_tmp = df_tmp[~df_tmp['gene'].isna()]
        gene_list = []
        for i in df_tmp.index:
            gene_list += df_tmp.loc[i, 'gene'].split(":")
        for gene in title_gene_list:
            df2.loc[area, gene] = gene_list.count(gene)
    return df2


def convert_in_hospital(input_df):
    df1 = copy.deepcopy(input_df)
    df1 = df1[df1['member_count'] == 2]
    columns=['hospital', 'reports', 'risk_couples', 'total_rate', 'auto_risk', 'auto_rate', 'x_couples', 'x_rate']
    hospital_list = list(set(df1['hospital'].tolist()))
    df3 = pd.DataFrame(columns=columns)
    hospital_count = 0
    for i in hospital_list:
        df_tmp = df1[df1['hospital'] == i]
        reports = df_tmp.shape[0]
        df3.loc[hospital_count, 'hospital'] = i
        df3.loc[hospital_count, 'reports'] = reports
        df3.loc[hospital_count, 'risk_couples'] = df_tmp[(df_tmp['x_gene_num'] > 0) |
                                                         (df_tmp['gene_num_at_risk'] > 0)].shape[0]
        df3.loc[hospital_count, 'total_rate'] = df3.loc[hospital_count, 'risk_couples'] / reports
        df3.loc[hospital_count, 'auto_risk'] = df_tmp[(df_tmp['gene_num_at_risk'] > 0)].shape[0]
        df3.loc[hospital_count, 'auto_rate'] = df3.loc[hospital_count, 'auto_risk'] / reports
        df3.loc[hospital_count, 'x_couples'] = df_tmp[(df_tmp['x_gene_num'] > 0)].shape[0]
        df3.loc[hospital_count, 'x_rate'] = df3.loc[hospital_count, 'x_couples'] / reports
        hospital_count += 1
    return df3


def data2plot_gene(input_df, cut_line=1/200, area=None):
    pre_df = copy.deepcopy(input_df)

    # 增加一行“total”  统计总人数
    if area:
        pre_df.loc['total'] = pre_df.loc[area]
    else:
        pre_df.loc['total'] = [sum(pre_df[t]) for t in pre_df.columns]
    male_counts = pre_df.loc['total', 'individuals_male']
    total_counts = pre_df.loc['total', 'individuals_total']
    female_counts = total_counts - male_counts
    glist = pre_df.columns.tolist()[5:]

    # 常隐和x连锁分别计算携带频率
    for gene in glist:
        if gene in Auto_list:
            pre_df.loc['total', gene] = pre_df.loc['total', gene] / total_counts
        elif gene in Xlink_list:
            pre_df.loc['total', gene] = pre_df.loc['total', gene] / female_counts
        else:
            raise ValueError

    # 去除人数相关的列，只保留基因，并以携带频率排序
    pre_df = pre_df[glist]
    pre_df.sort_values(by='total', inplace=True, axis=1)
    if cut_line:
        pre_df = pre_df.T[pre_df.T['total'] > cut_line].T
    return pre_df


def transform_area_gene_cf_matrix(input_df, cut_line=1/200):
    if cut_line:
        pre_df = fst.filter_by_cf(input_df, cut_line)
    else:
        pre_df = copy.deepcopy(input_df)
    glist = pre_df.columns.tolist()[5:]

    # 分别计算常隐和x连锁基因的携带频率
    for i in pre_df.index:
        for gene in glist:
            if gene in Auto_list:
                pre_df.loc[i, gene] = (pre_df.loc[i, gene]) / (pre_df.loc[i, 'individuals_total'])
            elif gene in Xlink_list:
                pre_df.loc[i, gene] = (pre_df.loc[i, gene]) / \
                                  (pre_df.loc[i, 'individuals_total'] - pre_df.loc[i, 'individuals_male'])
            else:
                raise ValueError
    return pre_df[glist]


def random_n_distance(n, input_df, area1='湖南'):
    df_sample = input_df
    df_sample_area1 = df_sample[df_sample.area == area1].sample(n, replace=False, axis=0)
    df = convert_in_areas(df_sample_area1)
    df.set_index('area', inplace=True)
    df = transform_area_gene_cf_matrix(df, 0)
    return np.linalg.norm(df)


def transform_merge_area(input_df, merge_rules):
    pre_df = copy.deepcopy(input_df)
    for i in merge_rules.keys():
        df_tmp = pre_df[[(t in merge_rules[i]) for t in pre_df.index]]      # 截取需合并地区数据
        pre_df.drop([t for t in merge_rules[i]], inplace=True)              # 删除原地区数据

        # 统计合并地区的携带者人数、总人数
        pre_df.loc[i, 'carriers_auto'] = sum(df_tmp['carriers_auto'].tolist())
        pre_df.loc[i, 'carriers_x'] = sum(df_tmp['carriers_x'].tolist())
        pre_df.loc[i, 'carriers_total'] = sum(df_tmp['carriers_total'].tolist())
        pre_df.loc[i, 'individuals_male'] = sum(df_tmp['individuals_male'].tolist())
        pre_df.loc[i, 'individuals_total'] = sum(df_tmp['individuals_total'].tolist())

        # 统计合并地区各个基因的检出次数
        l_cumulator = df_tmp.columns.tolist()
        del l_cumulator[0:5]
        for t in l_cumulator:
            pre_df.loc[i, t] = sum(df_tmp[t].tolist())
    return pre_df


def transform_data_for_stats(input_df, df_id):
    column = ['name', 'sex', 'id', 'fid', 'area', 'second_area', 'main_area'] + Auto_list + Xlink_list
    arr = np.zeros(input_df.shape)
    df_sample = pd.DataFrame(arr, columns=column)

    # 复制共同column到新df
    df_sample.name = input_df.name
    df_sample.sex = input_df.sex
    df_sample.fid = input_df.fid
    df_sample.area = input_df.area

    # 获取main_area和second_area
    main_area = [get_keys(area_counterparts2, i) for i in df_sample['area']]
    second_area = [get_keys(area_counterparts, i) for i in df_sample['area']]
    id_list = [get_sample_id(a, b, df_id) for a, b in zip(input_df.name, input_df.fid)]
    df_sample['main_area'] = main_area
    df_sample['second_area'] = second_area
    df_sample['id'] = id_list

    # 标注检出变异的基因
    for i in input_df.index:
        if input_df.loc[i, 'carrier_status']:
            glist = input_df.loc[i, 'gene'].split(':')
            for t in glist:
                df_sample.loc[i, t] = 1
    return df_sample


def province_ch_to_en_index(df, type_='index'):
    df_trans = pd.read_excel('D:\我的坚果云\ECS_1.6w_samples\province translation.xlsx', index_col='ch')
    if type_ == 'index':
        for i in df.index:
            df.rename(index={i: df_trans.loc[i, 'en']}, inplace=True)
    elif type_ == 'column':
        for i in df.columns:
            df.rename(columns={i: df_trans.loc[i, 'en']}, inplace=True)
    elif type_ == 'both':
        for i in df.columns:
            df.rename(columns={i: df_trans.loc[i, 'en']}, inplace=True)
        for i in df.index:
            df.rename(index={i: df_trans.loc[i, 'en']}, inplace=True)
    else:
        raise ValueError
