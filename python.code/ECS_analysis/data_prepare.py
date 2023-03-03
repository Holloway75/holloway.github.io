from addresss import *
import fst_calculation as fst


def get_keys(dic, val):
    for k, v in dic.items():
        if val in v:
            return k
    return val


def get_areas_from_id(name, fid, data):
    df1 = data[(data.fid == fid) & (data.name == name)]
    if not df1.shape[0]:
        return "unknown"
    elif df1.shape[0] > 1:
        raise ValueError
    else:
        return data.loc[df1.index[0], 'id']


def convert_in_areas(df):
    constant_column = ['area', 'carriers_auto', 'carriers_x', 'carriers_total',
                       'individuals_male', 'individuals_total']
    title_gene_list = Auto_list + Xlink_list
    column = constant_column + title_gene_list
    df2 = pd.DataFrame(columns=column)
    areas = list(set(df.area.tolist()))
    province_count = 0
    for x in areas:
        df2.loc[province_count, 'area'] = x
        df_x = df[df.area == x]

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
        for gene in title_gene_list:
            if pd.isna(df2.loc[province_count, gene]):
                df2.loc[province_count, gene] = 0

        province_count += 1
    return df2


def data2plot_gene(input_df, cut_line=1/200, area=None):
    if area is not None:
        pre_df = input_df.loc[area]
    else:
        pre_df = copy.deepcopy(input_df)
    if cut_line:
        pre_df = fst.filter_by_cf(pre_df, cut_line)

    # 增加一行“total”  统计总人数
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
    main_area = [get_keys(Area_counterparts2, i) for i in df_sample['area']]
    second_area = [get_keys(Area_counterparts, i) for i in df_sample['area']]
    id_list = [get_areas_from_id(a, b, df_id) for a, b in zip(input_df.name, input_df.fid)]
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


def transform_in_samples(df, data, sex_label):
    """
    将原始数据中所有样本性别相同，转录到新df2, 每行为一个个体
    :param df: 原始数据
    :param data: 包含name, fid与id(身份证前六位)
    :param sex_label: 原始数据中样本的性别，int，0 or 1
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

        auto_status, x_status, f8_inv_status, fmr1_status = 0, 0, 0, 0
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

        if sum([(x in Auto_list) for x in sample_gene_list]):
            auto_status = 1
        if sum([(x in Xlink_list) for x in sample_gene_list]):
            x_status = 1
        carrier_status = auto_status*2**3 + x_status*2**2 + f8_inv_status*2 + fmr1_status
        df2.loc[sample_counts, 'carrier_status'] = carrier_status
        df2.loc[sample_counts, 'gene'] = ":".join(sample_gene_list)
        df2.loc[sample_counts, 'var_id'] = ":".join(sample_var_list)
        lines += t
        sample_counts += 1
    return df2
