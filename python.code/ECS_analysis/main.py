from addresss import *
from plot import *
<<<<<<< HEAD
import matplotlib.pyplot as plt
=======
>>>>>>> f85d5d3f0f61e82a369e5026d6023a7710853ed9


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


<<<<<<< HEAD
def convert_in_areas(df, ):
    constant_column = ['area', 'carriers_auto', 'carriers_x', 'carriers_total',
                       'individuals_male', 'individuals_total']
    title_genelist = auto_list + xlink_list
=======
def convert_in_areas(df, alist, xlist):
    constant_column = ['area', 'carriers_auto', 'carriers_x', 'carriers_total',
                       'individuals_male', 'individuals_total']
    title_genelist = alist + xlist
>>>>>>> f85d5d3f0f61e82a369e5026d6023a7710853ed9
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


if __name__ == '__main__':
    df_area = pd.read_csv('area.combined.csv', index_col='area')
<<<<<<< HEAD
    df_merge = transform_merge_area(df_area, area_counterparts2)
    df_merge.drop('unknown', inplace=True)

    # plot_area2_fst_clustermap(df_merge)
    a = ['北方', '南方', '华南']
    fig ,ax = plt.subplots(1,3, figsize=(24,6))
    cut_line = 1/200
    for i in range(3):
        pre_df = data_prepare.data2plot_gene(df_merge, area=[a[i]])
        y = pre_df.loc['total']
        gene_num = len(y)
        x = np.arange(gene_num)

        ax[i].barh(x, y, height=0.7, color='#1f77b4', align='center', tick_label=pre_df.columns.tolist())
        ax[i].set_xscale("log")
        ax[i].spines['top'].set_color(None)
        ax[i].spines['right'].set_color(None)
        ax[i].set_title('Carrier frequency distribution of %d filtered genes' % gene_num, fontsize=14)

        plt.legend(loc=0, fontsize=12)
        ax[i].set_xlabel('Carrier frequency', fontsize=12)
        ax[i].set_ylabel('Gene', fontsize=12)
        if gene_num > 60:
            ax[i].set_yticks([])  # 基因数大于60，不显示基因名
    plt.show()

    df_merge = transform_merge_area(df_area, area_counterparts)
    df_merge.drop('unknown', inplace=True)

    plot_area2_fst_clustermap(df_merge)



