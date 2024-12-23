import os
from cProfile import label
from pickle import PROTO

import matplotlib.pyplot
from openpyxl.styles.builtins import title
from patsy.user_util import demo_data
from scipy.stats import alpha

import data_prepare
import statsmodels.api as sm
os.environ["OMP_NUM_THREADS"] = '1'
os.environ["KERAS_BACKEND"] = "jax"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':
    os.chdir('E:\我的坚果云')
    plt.rcParams['font.family'] = 'SimHei'
    plt.rcParams['font.size'] = 12

    df = pd.read_excel('wxy.xlsx', index_col=False)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes((0.15, 0.3, 0.8, 0.65))
    df = df[df.index < 9]

    sns.scatterplot(data=df, x='孕周', y='血小板（*10^9/L）')
    sns.lineplot(data=df, x='孕周', y='血小板（*10^9/L）')
    ax.hlines(y=70,xmin=0, xmax=35, colors='red', linestyles='--')
    ax.hlines(y=100, xmin=0, xmax=35, colors=sns.color_palette(n_colors=3)[1], linestyles='--')

    for i, row in df.iterrows():
        plt.text(row['孕周'], row['血小板（*10^9/L）'], f'{row["血小板（*10^9/L）"]}', ha='right', va='bottom')

    # ax.set_xticks(np.arange(8,15))
    # ax.set_xticklabels(df['label'], rotation=45, ha='right', va='top')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()
