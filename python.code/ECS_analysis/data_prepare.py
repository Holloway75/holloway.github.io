from addresss import *
import fst_calculation as fst


def convert_in_areas(df):
    constant_column = ['area', 'carriers_auto', 'carriers_x', 'carriers_total',
                       'individuals_male', 'individuals_total']
    title_genelist = auto_list + xlink_list
    column = constant_column + title_genelist
    df2 = pd.DataFrame(columns=column)
    areas = list(set(df.area.tolist()))
    province_count = 0
    for x in areas:
        df2.loc[province_count, 'area'] = x
        df_x = df[df.area==x]

        # 计数该地区中 常染色体病携带者 x携带者 总携带者人数 男性人数 总人数
        status_list = df_x.carrier_status.tolist()
        df2.loc[province_count, 'carriers_auto'] = len([sam for sam in status_list if sam >= 8])
        df2.loc[province_count, 'carriers_x'] = len([sam for sam in status_list if (sam % 8) >= 4])
        df2.loc[province_count, 'carriers_total'] = len([sam for sam in status_list if sam > 0])
        df2.loc[province_count, 'individuals_male'] = sum(df_x.sex.tolist())
        df2.loc[province_count, 'individuals_total'] = len(status_list)

        # 计数该地区中每个基因的出现次数
        raw_str_list = df_x.gene.tolist()
        str_list = [gene for gene in raw_str_list if isinstance(gene, str)]         # 去除nan
        gene_list = []
        for strs in str_list:
            for single_gene in strs.split(':'):         # 拆分字符串形式的基因“gene1:gene2：gene3”
                gene_list.append(single_gene)
        gene_list_rmdup = list(set(gene_list))
        for gene in gene_list_rmdup:
            df2.loc[province_count, gene] = gene_list.count(gene)

        # 该省未检出的基因填0
        for gene in title_genelist:
            if pd.isna(df2.loc[province_count, gene]):
                df2.loc[province_count, gene] = 0

        province_count += 1
    return df2


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


def transform_area_gene_cf_matrix(input_df, cut_line=1/200):
    if cut_line:
        pre_df = fst.filter_by_cf(input_df, cut_line)
    else:
        pre_df = copy.deepcopy(input_df)
    glist = pre_df.columns.tolist()[5:]

    # 分别计算常隐和x连锁基因的携带频率
    for i in pre_df.index:
        for gene in glist:
            if gene in auto_list:
                pre_df.loc[i, gene] = (pre_df.loc[i, gene]) / (pre_df.loc[i, 'individuals_total'])
            elif gene in xlink_list:
                pre_df.loc[i, gene] = (pre_df.loc[i, gene]) / \
                                  (pre_df.loc[i, 'individuals_total'] - pre_df.loc[i, 'individuals_male'])
            else:
                raise ValueError
    return pre_df[glist]


def random_n_distance(n, input_df, area1='湖南'):
    df_sample = input_df
    df_sample_area1 = df_sample[df_sample.area==area1].sample(n, replace=False, axis=0)
    df = convert_in_areas(df_sample_area1)
    df.set_index('area', inplace=True)
    df = transform_area_gene_cf_matrix(df, 0)
    return np.linalg.norm(df)



