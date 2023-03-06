import matplotlib.pyplot as plt
import pandas as pd

import stats
from addresss import *
from plot import *


if __name__ == '__main__':
    df = pd.read_csv('sample.stats.csv')
    df_stats = pd.DataFrame(columns=['gene', 'case', 'control', 'odds_ratio', 'ci_lower', 'ci_upper', 'rank'])
    for i in Gene_over_200:
        index = Gene_over_200.index(i)
        r1 = stats.calculate_odds_ratio(df, control_label='北方', case_label='南方', label_column='main_area', gene=i)
        r2 = stats.calculate_odds_ratio(df, control_label='北方', case_label='华南', label_column='main_area', gene=i)
        df_stats.loc[2*index] = [i, r1.case, r1.control, r1.value, r1.interval[0], r1.interval[1], (r1.value+r2.value)]
        df_stats.loc[2*index+1] = [i, r2.case, r2.control, r2.value, r2.interval[0], r2.interval[1], (r1.value+r2.value)]

    df_stats.sort_values(by='rank', inplace=True, ignore_index=True)

    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    df_1 = df_stats[df_stats.case == '南方'].reset_index(drop=True)
    df_2 = df_stats[df_stats.case == '华南'].reset_index(drop=True)

    for i in range(4):
        df_tmp1 = df_1.loc[7*i:(7*i+6)]
        df_tmp2 = df_2.loc[7 * i:(7 * i + 6)]
        ax[i//3][i % 3].errorbar(x=df_tmp1['odds_ratio'],
                                 y=np.arange(0, 14, 2),
                                 xerr=[(df_tmp1.odds_ratio - df_tmp1['ci_lower']),
                                       (df_tmp1.ci_upper - df_tmp1.odds_ratio)],
                                 fmt='o', mfc='red', mec='red', capsize=3, markersize=5, ecolor='tab:orange')
        ax[i // 3][i % 3].errorbar(x=df_tmp2['odds_ratio'],
                                   y=np.arange(-0.5, 12, 2),
                                   xerr=[(df_tmp2.odds_ratio - df_tmp2['ci_lower']),
                                         (df_tmp2.ci_upper - df_tmp2.odds_ratio)],

                                   fmt='o', mfc='red', mec='red', capsize=3, markersize=5, ecolor='tab:green')
        ax[i//3][i % 3].spines['top'].set_color(None)
        ax[i//3][i % 3].spines['right'].set_color(None)
        ax[i//3][i % 3].spines['left'].set_visible(False)
        ax[i//3][i % 3].spines['bottom'].set_position(('data', -1))
        ax[i//3][i % 3].tick_params(left=False)
        ax[i//3][i % 3].vlines(1, -1, 13, color='black', lw=1)
        ax[i // 3][i % 3].set_yticks(np.arange(-0.5, 12, 2)+0.25, df_tmp1.gene.tolist())

    i = 4
    df_tmp1 = df_1.loc[7 * i:(7 * i + 5)]
    df_tmp2 = df_2.loc[7 * i:(7 * i + 5)]
    ax[i // 3][i % 3].errorbar(x=df_tmp1['odds_ratio'],
                               y=np.arange(0, 12, 2),
                               xerr=[(df_tmp1.odds_ratio - df_tmp1['ci_lower']),
                                     (df_tmp1.ci_upper - df_tmp1.odds_ratio)],
                               fmt='o', mfc='red', mec='red', capsize=3, markersize=5, ecolor='tab:orange')
    ax[i // 3][i % 3].errorbar(x=df_tmp2['odds_ratio'],
                               y=np.arange(-0.5, 10, 2),
                               xerr=[(df_tmp2.odds_ratio - df_tmp2['ci_lower']),
                                     (df_tmp2.ci_upper - df_tmp2.odds_ratio)],
                               fmt='o', mfc='red', mec='red', capsize=3, markersize=5, ecolor='tab:green')
    ax[i // 3][i % 3].spines['top'].set_color(None)
    ax[i // 3][i % 3].spines['right'].set_color(None)
    ax[i // 3][i % 3].spines['left'].set_visible(False)
    ax[i // 3][i % 3].spines['bottom'].set_position(('data', -1))
    ax[i // 3][i % 3].tick_params(left=False)
    ax[i // 3][i % 3].vlines(1, -1, 11, color='black', lw=1)
    ax[i // 3][i % 3].set_yticks(np.arange(-0.5, 11, 2) + 0.25, df_tmp1.gene.tolist())

    i = 5
    df_tmp1 = df_stats[(df_stats.gene == 'G6PD') & (df_stats.case == '南方')]
    df_tmp2 = df_stats[(df_stats.gene == 'G6PD') & (df_stats.case == '华南')]
    ax[i // 3][i % 3].errorbar(x=df_tmp1['odds_ratio'],
                               y=[4],
                               xerr=[(df_tmp1.odds_ratio - df_tmp1['ci_lower']),
                                     (df_tmp1.ci_upper - df_tmp1.odds_ratio)],
                               fmt='o', mfc='red', mec='red', capsize=3, markersize=5, ecolor='tab:orange')
    ax[i // 3][i % 3].errorbar(x=df_tmp2['odds_ratio'],
                               y=[2],
                               xerr=[(df_tmp2.odds_ratio - df_tmp2['ci_lower']),
                                     (df_tmp2.ci_upper - df_tmp2.odds_ratio)],
                               fmt='o', mfc='red', mec='red', capsize=3, markersize=5, ecolor='tab:green')
    ax[i // 3][i % 3].spines['top'].set_color(None)
    ax[i // 3][i % 3].spines['right'].set_color(None)
    ax[i // 3][i % 3].spines['left'].set_visible(False)
    ax[i // 3][i % 3].spines['bottom'].set_position(('data', -1))
    ax[i // 3][i % 3].tick_params(left=False)
    ax[i // 3][i % 3].vlines(1, -1, 11, color='black', lw=1)
    ax[i // 3][i % 3].set_yticks([3], ['G6PD'])


    plt.show()


