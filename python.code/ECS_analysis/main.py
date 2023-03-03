import pandas as pd

import stats
from addresss import *
from plot import *


if __name__ == '__main__':
    df = pd.read_csv('sample.stats.csv')
    df_stats = pd.DataFrame(columns=['gene', 'case_label', 'odds_ratio', 'ci'])
    i = 'GJB2'
    r1 = stats.calculate_odds_ratio(df, control_label='北方', case_label='南方', label_column='main_area', gene=i)
    r2 = stats.calculate_odds_ratio(df, control_label='北方', case_label='华南', label_column='main_area', gene=i)
    df_stats.loc[0] = [i, r1.case, r1.value, r1.interval]
    df_stats.loc[1] = [i, r2.case, r2.value, r2.interval]

    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(x=[2, 4], y=df_stats.odds_ratio, fmt='o')

    plt.show()





