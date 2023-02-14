import matplotlib.pyplot as plt
import seaborn as sns
import fst_calculation as fst
from addresss import *
from pyecharts import options as opts
from pyecharts.charts import Map
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from fractions import Fraction
from sklearn.metrics import silhouette_score, calinski_harabasz_score



def std_pca(data, n_components=2):
    st = StandardScaler()
    data_std = st.fit_transform(data)
    pca = PCA(n_components=n_components)
    data_pca = pca.fit(data_std).transform(data_std)
    return data_pca


def plot_china_map(input_df):
    data = input_df[['area', 'individuals']].values.tolist()
    for i in data:
        i[0] = province_name_simple_to_full(i[0])

    c = Map()
    c.add(series_name='样本数量', data_pair=data, maptype="china")
    c.set_global_opts(title_opts=opts.TitleOpts(title="样本收集数量"), visualmap_opts=opts.VisualMapOpts(max_=1000))
    c.set_series_opts(label_opts=opts.LabelOpts(is_show=False))
    c.render("map_base.html")


def plot_area_individual(input_df):
    pre_df = copy.deepcopy(input_df)
    pre_df.drop('unknown', inplace=True)
    pre_df.sort_values('individuals_total', inplace=True, ascending=False)

    fig, ax = plt.subplots(figsize=(8, 6))
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

    x = np.arange(pre_df.shape[0])
    ax.bar(x, pre_df.carriers_total, width=0.4, edgecolor="white", color='#E0884B', align='edge', label='常染色体携带者')
    ax.bar(x+0.4, pre_df.individuals_total, width=0.4, edgecolor="white", tick_label=pre_df.index,
           color='#1f77b4', align='edge', label='接受筛查人数')

    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    plt.axhline(100, color='r', label='样本量=100')
    ax.set_xlabel('地区',  fontsize=16)
    ax.set_ylabel('人数', fontsize=16)
    ax.set_title('各地区样本量及常染色体携带者人数', fontsize=20)
    plt.legend(loc="upper right", fontsize=16)
    plt.show()


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
        l = df_tmp.columns.tolist()
        del l[0:5]
        for t in l:
            pre_df.loc[i, t] = sum(df_tmp[t].tolist())
    return pre_df


def plot_gene(input_df, cut_line=1/500, picture=True):
    pre_df = copy.deepcopy(input_df)

    # 增加一行'total' 并统计总人数
    pre_df.loc['total'] = [sum(pre_df[t]) for t in pre_df.columns]
    male_counts = pre_df.loc['total', 'individuals_male']
    total_counts = pre_df.loc['total', 'individuals_total']
    female_counts = total_counts - male_counts

    # 常隐和x连锁分别计算携带频率
    for gene in pre_df.columns.tolist()[5:]:
        if gene in auto_list:
            pre_df.loc['total', gene] = pre_df.loc['total', gene] / total_counts
        elif gene in xlink_list:
            pre_df.loc['total', gene] = pre_df.loc['total', gene] / female_counts
        else:
            raise ValueError

    # 去除人数相关的列，只保留基因，并以携带频率排序
    pre_df.drop(['carriers_auto', 'carriers_x', 'carriers_total', 'individuals_male', 'individuals_total'],
                inplace=True, axis=1)
    pre_df.sort_values(by='total', inplace=True, axis=1)

    # 绘图
    xtickslabel = [t for t in pre_df.columns if pre_df.loc['total', t] > cut_line]
    if picture:
        fig, ax = plt.subplots(figsize=(8, 6))
        y = [t for t in pre_df.loc['total'] if t > cut_line]
        gene_num = len(y)
        x = np.arange(gene_num)

        ax.barh(x, y, height=0.7, color='#1f77b4', align='center', tick_label=xtickslabel)
        ax.set_xscale("log")
        ax.spines['top'].set_color(None)
        ax.spines['right'].set_color(None)
        if cut_line:
            plt.axvline(cut_line, color='r', label='Carrier frequency=%s' % str(Fraction(1, int(1/cut_line))))
            ax.set_title('Carrier frequency distribution of %d filtered genes' % gene_num, fontsize=14)

        else:
            plt.axvline(1/500, color='r', label='Carrier frequency=1/500')
            ax.set_title('Carrier frequency distribution of %d genes' % gene_num, fontsize=14)
        plt.legend(loc=0, fontsize=12)
        ax.set_xlabel('Carrier frequency', fontsize=12)
        ax.set_ylabel('Gene', fontsize=12)
        if gene_num > 60:
            ax.set_yticks([])       # 基因数大于60，不显示基因名
        plt.show()
    return xtickslabel          # 返回保留的基因列表


def transform_area_gene_cf_matrix(input_df, cut_line=1/200):
    pre_df = copy.deepcopy(input_df)
    glist = plot_gene(pre_df, cut_line, picture=False)

    # 分别计算常隐和x连锁基因的携带频率
    for i in pre_df.index:
        for gene in auto_list:
            pre_df.loc[i, gene] = (pre_df.loc[i, gene]+1) / (pre_df.loc[i, 'individuals_total']+2)
        for gene in xlink_list:
            pre_df.loc[i, gene] = (pre_df.loc[i, gene]+1) / \
                                  (pre_df.loc[i, 'individuals_total'] - pre_df.loc[i, 'individuals_male']+2)
    return pre_df[glist]


def plot_area_pca(input_df, cut_line=1/200, std=True):
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    pre_df = transform_area_gene_cf_matrix(input_df=input_df, cut_line=cut_line)

    fig, ax = plt.subplots(figsize=(8, 6))
    x = np.array(pre_df)
    if std:
        ax.set_title('各地区单基因病携带频率主成分分析--标准化', fontsize=14)
        x_r = std_pca(x, 2)
    else:
        ax.set_title('各地区单基因病携带频率主成分分析--未准化', fontsize=14)
        pca = PCA(n_components=2)
        x_r = pca.fit(x).transform(x)
    ax.scatter(x_r[:, 0], x_r[:, 1], alpha=0.9, lw=2)

    # 标注地域名称
    count = 0
    for text in pre_df.index:
        if text == '京津':
            plt.annotate(text, xy=(x_r[count, 0], x_r[count, 1]), xytext=(x_r[count, 0], x_r[count, 1] + 0.001),
                         fontsize=12)
        elif text == '安徽':
            plt.annotate(text, xy=(x_r[count, 0], x_r[count, 1]), xytext=(x_r[count, 0] - 0.001, x_r[count, 1] - 0.001),
                         fontsize=12)
        else:
            plt.annotate(text, xy=(x_r[count, 0], x_r[count, 1]), xytext=(x_r[count, 0] + 0.001, x_r[count, 1] - 0.001),
                         fontsize=12)
        count += 1

    ax.set_xlabel('PC1', fontsize=12)
    ax.set_ylabel('PC2', fontsize=12)
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    plt.show()


def plot_kmeans_pca(input_df, k=2, picture=True, random_state=150):
    x = np.array(input_df)
    y_pred = KMeans(n_clusters=k, random_state=random_state).fit_predict(x)

    if picture:
        x_r = std_pca(x, 2)
        fig, ax = plt.subplots(figsize=(8, 6))
        for i in range(k):
            ax.scatter(x_r[y_pred == i, 0], x_r[y_pred == i, 1], alpha=0.8, lw=2)
        for text, count in zip(input_df.index, np.arange(input_df.shape[0])):
            if text == '京津':
                plt.annotate(text, xy=(x_r[count, 0], x_r[count, 1]), xytext=(x_r[count, 0], x_r[count, 1] + 0.001),
                             fontsize=12)
            elif text == '安徽':
                plt.annotate(text, xy=(x_r[count, 0], x_r[count, 1]), xytext=(x_r[count, 0] - 0.001, x_r[count, 1] - 0.001),
                             fontsize=12)
            else:
                plt.annotate(text, xy=(x_r[count, 0], x_r[count, 1]), xytext=(x_r[count, 0] + 0.001, x_r[count, 1] - 0.001),
                             fontsize=12)
        plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
        plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
        ax.set_xlabel('PC1', fontsize=12)
        ax.set_ylabel('PC2', fontsize=12)
        ax.set_title('各地区单基因病携带频率主成分分析(K=%d)' % k, fontsize=14)
        ax.spines['top'].set_color(None)
        ax.spines['right'].set_color(None)
        plt.show()
    return [silhouette_score(x, y_pred), calinski_harabasz_score(x, y_pred)]


def kmeans_evaluations(input_df):
    sil, cal = [], []
    for i in range(2, 11):
        evaluations_index_list = plot_kmeans_pca(input_df, k=i, picture=False)
        sil.append(evaluations_index_list[0])
        cal.append(evaluations_index_list[1])
    fig, ax1 = plt.subplots(figsize=(8, 6))
    x = [str(i) for i in range(2, 11)]
    ax1.plot(x, sil,  'bo-', label='silhouette score')
    ax2 = ax1.twinx()
    ax2.plot(x, cal, 'ro-', label='calinski harabasz score')
    fig.legend(labels = ('Silhouette','Calinski Harabasz'), loc=(0.65, 0.75))
    ax1.set_xlabel('Groups', fontsize=12)
    ax1.set_ylabel('Silhouette', fontsize=12)
    ax2.set_ylabel('Calinski Harabasz', fontsize=12)
    ax1.set_title('KMeans Evaluation', fontsize=14)
    plt.show()


def plot_area_gene_heatmap(input_df):
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    pre_df = copy.deepcopy(input_df)
    pre_df.sort_values(by='GJB2', axis=0, inplace=True)
    area_heatmap_list = pre_df.index.tolist()
    gene_heatmap_list = pre_df.columns.tolist()
    x = np.array(pre_df)
    x = np.log(x) / np.log(10)
    fig, ax = plt.subplots(figsize=(16, 9))
    im = ax.imshow(x)
    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(gene_heatmap_list)), labels=gene_heatmap_list)
    ax.set_yticks(np.arange(len(area_heatmap_list)), labels=area_heatmap_list)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax.set_title("各地区常见单基因病携带频率分布差异")
    fig.tight_layout()
    plt.show()


def plot_area_area_heatmap(input_df):
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    input_df = transform_area_gene_cf_matrix(input_df, cut_line=1/200)
    area_heatmap_list = input_df.index.tolist()

    data = np.array(input_df)
    st = StandardScaler()
    data_std = st.fit_transform(data)
    # 计算各个地区之间距离
    pre_df = pd.DataFrame()

    for arr in area_heatmap_list:
        for col in area_heatmap_list:
            arr_num = area_heatmap_list.index(arr)
            col_num = area_heatmap_list.index(col)
            pre_df.loc[arr, col] = np.linalg.norm([b-a for b,a in zip(data_std[arr_num], data_std[col_num])])

    pre_df.index = pre_df.index.astype('category').set_categories(area_sort_list, ordered=True)
    pre_df.sort_index(inplace=True)
    pre_df = pre_df[area_sort_list]

    area_heatmap_list = pre_df.index.tolist()

    x = np.array(pre_df)
    distance_max = np.max(x)
    x = 1 - x/distance_max

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(x)
    # Show all ticks and label them with the respective list entries
    ax.set_xticks(np.arange(len(area_heatmap_list)), labels=area_heatmap_list)
    ax.set_yticks(np.arange(len(area_heatmap_list)), labels=area_heatmap_list)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax.set_title("各地区单基因病携带谱差异性分析")
    fig.tight_layout()
    plt.show()


def plot_area2_fst_heatmap(input_df):
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    pre_df2 = fst.data_prepare_for_heatmap(input_df)

    mask = np.triu(np.ones_like(pre_df2, dtype=bool), 1)    # 遮盖上三角
    fig, ax = plt.subplots(figsize=(16, 9))
    sns.heatmap(pre_df2, mask=mask, cmap='YlGnBu_r', robust=True, annot=True,
                annot_kws={'size': 9, 'weight':'bold'},
                fmt='.4f', square=True, linewidths=.5, cbar_kws={"shrink": .5})

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    ax.set_title("各地区携带频率差异性-Fst")
    fig.tight_layout()
    plt.show()


def plot_area2_fst_clustermap(input_df):
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    pre_df2 = fst.data_prepare_for_heatmap(input_df)


    g = sns.clustermap(
        pre_df2,
        figsize=(8, 6),
        row_cluster=False,
        dendrogram_ratio=(0.3, 0.2),
        cbar_pos = (0.1, .2, .03, .4),
        cmap='YlGnBu_r', robust=True
    )

    plt.show()


