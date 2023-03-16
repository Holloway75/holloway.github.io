import numpy as np
import pandas as pd
import addresss
import data_prepare
import plot


if __name__ == '__main__':
    df = pd.read_csv('area.merged.csv', index_col='area')
    df.drop(index='unknown', inplace=True)
    plot.plot_area2_fst_clustermap(df)