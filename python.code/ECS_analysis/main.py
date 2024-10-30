import os
from cProfile import label
from pickle import PROTO

import matplotlib.pyplot
from openpyxl.styles.builtins import title
from scipy.stats import alpha

import data_prepare
import statsmodels.api as sm
os.environ["OMP_NUM_THREADS"] = '1'
os.environ["KERAS_BACKEND"] = "jax"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from addresss import *


def plot_model(ax_: matplotlib.pyplot.axes):
    os.chdir('D:\我的坚果云\投稿-卫生经济学')
    df = pd.read_excel('data.xlsx', index_col=0, sheet_name='gene_count')
    df.sort_values(by='lng', inplace=True)

    x_ = df[['lng', 'c', 's']]
    y_ = df['y']

    x_ = sm.add_constant(x_)
    model = sm.OLS(y_, x_).fit()
    x_simu = np.ones((100, 4))
    x_simu[:, 1] = np.log(np.linspace(3, 1400, 100))
    x_seq = np.ones((100, 4))
    x_seq[:, 1] = np.log(np.linspace(3, 1400, 100))
    x_seq[:, 2] = 0
    y_simu = model.predict(x_simu)
    y_seq = model.predict(x_seq)

    sns.scatterplot(data=df, x='g', y='y', hue='c', s=15, hue_order=[1, 0])
    sns.lineplot(
        x=np.linspace(3, 1400, 100),
        y=y_simu,
        color=sns.color_palette(n_colors=2)[0],
        linestyle='--',
        label='Simultaneously',
        alpha=0.7
    )

    sns.lineplot(
        x=np.linspace(3, 1400, 100),
        y=y_seq,
        color=sns.color_palette(n_colors=2)[1],
        linestyle='--',
        label='Sequentially',
        alpha=0.7
    )

    ax_.set_ylabel('Unit Price, CHY')
    ax_.set_xlabel('Gene Counts in Panel')
    ax_.spines['top'].set_visible(False)
    ax_.spines['right'].set_visible(False)
    ax_.set_title('a', x=0.5, y=-0.3)

    # 提取图例句柄和标签
    handles, labels = ax_.get_legend_handles_labels()
    labels[0:2] = ['Simultaneously', 'Sequentially']
    ax.legend(handles, labels, title='Screening Method', markerscale=0.5)


def plot_arc(ax_: matplotlib.pyplot.axes):
    os.chdir('D:\我的坚果云\投稿-卫生经济学')
    df = pd.read_excel('data.xlsx', index_col=0, sheet_name='Sheet1')
    sns.lineplot(data=df, x='g', y='accu-arcr')
    l1 = [25, 0.0122767040317023]
    l2 = [48, 0.0135600046915091]
    ax_.vlines(l1[0], ymin=0, ymax=l1[1], label=f'1/100, {l1[0]} genes', colors=sns.color_palette(n_colors=4)[-2],
               linestyles='--')
    ax_.hlines(l1[1], xmin=0, xmax=l1[0], colors=sns.color_palette(n_colors=4)[-2], linestyles='--')
    ax_.vlines(l2[0], ymin=0, ymax=l2[1], label=f'1/200, {l2[0]} genes', colors=sns.color_palette(n_colors=4)[-1],
               linestyles='--')
    ax_.hlines(l2[1], xmin=0, xmax=l2[0], colors=sns.color_palette(n_colors=4)[-1], linestyles='--')
    ax_.set_ylim(0.002, 0.0145)
    ax_.set_xlim(0, 200)
    ax_.legend(title='Corresponding \n CF Threshold')
    ax_.spines['top'].set_visible(False)
    ax_.spines['right'].set_visible(False)
    ax_.set_xlabel('Gene Counts in Panel')
    ax_.set_ylabel('Combined At-Risk Couple Rate')
    ax_.set_title('b', x=0.5, y=-0.3)


def plot_curve(ax_: matplotlib.pyplot.axes):
    os.chdir('D:\我的坚果云\投稿-卫生经济学')
    df = pd.read_excel('data.xlsx', index_col=0, sheet_name='Sheet1')
    sns.lineplot(data=df, x='g', y='cost-effect-simu', label='Simultaneously')
    sns.lineplot(data=df, x='g', y='cost-effect-seq', label='Sequentially')

    ax_.set_xlabel('Gene Counts in Panel')
    ax_.set_title('c', x=0.5, y=-0.3)
    ax_.set_ylabel('Cost / Effectiveness, CHY')
    ax_.set_xlim(4, 200)
    ax_.set_ylim(0, 330000)
    l1 = [25, 205066.3953302]
    l2 = [48, 232370.740286589]
    ax_.vlines(l1[0], ymin=0, ymax=l1[1], label=f'1/100, {l1[0]} genes', colors=sns.color_palette(n_colors=4)[-2],
               linestyles='--')
    ax_.hlines(l1[1], xmin=0, xmax=l1[0], colors=sns.color_palette(n_colors=4)[-2], linestyles='--')
    ax_.vlines(l2[0], ymin=0, ymax=l2[1], label=f'1/200, {l2[0]} genes', colors=sns.color_palette(n_colors=4)[-1],
               linestyles='--')
    ax_.hlines(l2[1], xmin=0, xmax=l2[0], colors=sns.color_palette(n_colors=4)[-1], linestyles='--')
    ax_.legend(loc='lower right')
    ax_.spines['top'].set_visible(False)
    ax_.spines['right'].set_visible(False)


if __name__ == '__main__':
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 8
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_axes((0.12, 0.61, 0.35, 0.35))
    plot_model(ax)
    ax2 = fig.add_axes((0.60, 0.61, 0.35, 0.35))
    plot_arc(ax2)
    ax3 = fig.add_axes((0.15, 0.12, 0.75, 0.35))
    plot_curve(ax3)
    fig.savefig(r'D:\我的坚果云\投稿-卫生经济学\tmp.png', dpi=600)
    plt.show()

