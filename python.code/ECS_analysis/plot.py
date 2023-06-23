import matplotlib.pyplot as plt
import seaborn as sns
import fst_calculation as fst
import data_prepare
from pyecharts import options as opts
from pyecharts.charts import Map
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from fractions import Fraction
from sklearn.metrics import silhouette_score, calinski_harabasz_score
import copy
import pandas as pd
import numpy as np
from addresss import *


def plot_china_map(input_df):
    data = input_df[['area', 'individuals_total']].values.tolist()
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


def plot_gene(input_df, cut_line=1/200, area=None):
    pre_df = data_prepare.data2plot_gene(input_df, cut_line, area)
    fig, ax = plt.subplots(figsize=(8, 6))
    y = pre_df.loc['total']
    gene_num = len(y)
    x = np.arange(gene_num)

    ax.barh(x, y, height=0.7, color='#1f77b4', align='center', tick_label=pre_df.columns.tolist())
    ax.set_xscale("log")
    ax.spines['top'].set_color(None)
    ax.spines['right'].set_color(None)
    if cut_line:
        plt.axvline(cut_line, color='r', label='Carrier frequency=%s' % str(Fraction(1, int(1 / cut_line))))
        ax.set_title('Carrier frequency distribution of %d filtered genes' % gene_num, fontsize=14)

    else:
        plt.axvline(1 / 200, color='r', label='Carrier frequency=1/200')
        ax.set_title('Carrier frequency distribution of %d genes' % gene_num, fontsize=14)
    plt.legend(loc=0, fontsize=12)
    ax.set_xlabel('Carrier frequency', fontsize=12)
    ax.set_ylabel('Gene', fontsize=12)
    if gene_num > 60:
        ax.set_yticks([])  # 基因数大于60，不显示基因名
    plt.show()


def plot_kmeans_pca(input_df: pd.DataFrame, k=3, cut_line=1/200, picture=True):
    pre_df = data_prepare.transform_area_gene_cf_matrix(input_df, cut_line)
    x = np.array(pre_df)
    y_pred = KMeans(n_clusters=k, random_state=3).fit_predict(x)

    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

    if picture:
        pca = PCA(n_components=2)
        x_r = pca.fit_transform(x)
        df = pd.DataFrame(index=pre_df.index, columns=['PC1', 'PC2', 'area'], data=np.hstack((x_r, y_pred.reshape(-1,1))))
        df.sort_values(by='PC1', inplace=True)
        current_palette = sns.color_palette()

        north_name = "far south: " + ",".join(df[df['area'] == 1].index)
        south_name = "north: " + ",".join(df[df['area'] == 0].index)
        far_sor = 'south: ' + ",".join(df[df['area'] == 2].index)
        df['area'].replace({1:north_name, 0:south_name, 2:far_sor}, inplace=True)
        plt.figure(figsize=(16, 9))
        sns.scatterplot(df, x='PC1', y='PC2', hue='area', palette=current_palette)

        for i in df.index:
            plt.annotate(i, xy=(df.loc[i, 'PC1'], df.loc[i, 'PC2']), fontsize=12)
        plt.show()
    else:
        return [silhouette_score(x, y_pred), calinski_harabasz_score(x, y_pred)]


def kmeans_evaluations(input_df, cut_line=1/200):
    sil, cal = [], []
    for i in range(2, 11):
        evaluations_index_list = plot_kmeans_pca(input_df, k=i, cut_line=cut_line, picture=False)
        sil.append(evaluations_index_list[0])
        cal.append(evaluations_index_list[1])
    fig, ax1 = plt.subplots(figsize=(8, 6))
    x = [str(i) for i in range(2, 11)]
    ax1.plot(x, sil,  'bo-', label='silhouette score')
    ax2 = ax1.twinx()
    ax2.plot(x, cal, 'ro-', label='calinski harabasz score')
    fig.legend(labels=('Silhouette','Calinski Harabasz'), loc=(0.65, 0.75))
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


def plot_area_area_heatmap(df):
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    input_df = data_prepare.transform_area_gene_cf_matrix(df, cut_line=1/200)
    area_heatmap_list = input_df.index.tolist()

    # 计算各个地区之间距离
    a = np.array(input_df)
    b = np.ones((a.shape[1], a.shape[0]))
    result = np.sqrt((np.dot(a**2, b) + np.dot(a**2, b).T - 2 * np.dot(a, a.T)))
    pre_df = pd.DataFrame(index=area_heatmap_list, columns=area_heatmap_list, data=result)

    pre_df.sort_values(by='桂琼', inplace=True, ascending=False)
    pre_df.sort_values(by='桂琼', axis=1, inplace=True, ascending=False)

    area_heatmap_list = pre_df.index.tolist()

    x = np.array(pre_df)
    distance_max = np.max(x)
    x = 1 - x/distance_max

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(x, cmap='coolwarm')
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
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(pre_df2, mask=mask, cmap='coolwarm_r', robust=True, annot=True,
                annot_kws={'size': 9, 'weight':'bold'},
                fmt='.3f', square=True, linewidths=.5, cbar_kws={"shrink": .5})

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    ax.set_title("各地区携带频率差异性-Fst")
    fig.tight_layout()
    plt.show()


def plot_area2_fst_clustermap(input_df):
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    pre_df2 = fst.data_prepare_for_heatmap(input_df)
    g = sns.clustermap(pre_df2, figsize=(8, 6), row_cluster=False, dendrogram_ratio=(0.2, 0.2), method='ward',
                   cbar_pos=(0.1, .2, .03, .4), cmap='coolwarm_r', robust=True)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)

    plt.show()

