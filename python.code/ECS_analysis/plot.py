import matplotlib.pyplot as plt
import seaborn as sns
import fst_calculation as fst
import data_prepare
from pyecharts import options as opts
from pyecharts.charts import Map
from sklearn.decomposition import PCA
from fractions import Fraction
import copy
import pandas as pd
import numpy as np
from adjustText import adjust_text
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
    return pre_df


def plot_area_pca(input_df: pd.DataFrame, cut_line=1/200):
    pre_df = data_prepare.transform_area_gene_cf_matrix(input_df, cut_line)
    x = np.array(pre_df)

    pca = PCA(n_components=2)
    x_r = pca.fit_transform(x)
    df = pd.DataFrame(index=pre_df.index, columns=['PC1', 'PC2'], data=x_r)

    for i in df.index:
        if i in area_counterparts.keys():
            df.loc[i, 'Main Area'] = data_prepare.get_keys(area_counterparts2, area_counterparts[i][0])
        else:
            df.loc[i, 'Main Area'] = data_prepare.get_keys(area_counterparts2, i)
    data_prepare.province_ch_to_en_index(df)
    df.sort_values(by='PC1', inplace=True)
    current_palette = sns.color_palette()

    fig = plt.figure(figsize=(8, 6))
    sns.scatterplot(df, x='PC1', y='PC2', hue='Main Area', palette=current_palette)
    new_texts = [plt.text(x_, y_, text, fontsize=8) for x_, y_, text in zip(df['PC1'], df['PC2'], df.index)]
    adjust_text(new_texts, arrowprops=dict(arrowstyle='-', color='grey', lw=1))
    # plt.legend()
    return fig


def plot_area2_fst_heatmap(input_df, cut_line=1/200):
    pre_df2 = fst.data_prepare_for_heatmap(input_df, cut_line)
    data_prepare.province_ch_to_en_index(pre_df2, type_='both')
    mask = np.triu(np.ones_like(pre_df2, dtype=bool), 1)    # 遮盖上三角
    fig, ax = plt.subplots(figsize=(12, 9))
    sns.heatmap(pre_df2, mask=mask, cmap='coolwarm_r', robust=True, annot=True,
                annot_kws={'size': 6, 'weight':'bold'},
                fmt='.3f', square=True, linewidths=.5, cbar_kws={"shrink": .5})

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    fig.tight_layout()
    return fig


def plot_area2_fst_clustermap(input_df):
    pre_df2 = fst.data_prepare_for_heatmap(input_df)
    data_prepare.province_ch_to_en_index(pre_df2, type_='both')
    g = sns.clustermap(pre_df2, figsize=(8, 6), row_cluster=False, dendrogram_ratio=(0.2, 0.2), method='ward',
                   cbar_pos=(0.1, .2, .03, .4), cmap='coolwarm_r', robust=True)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)

    plt.show()

