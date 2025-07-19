import os
from tqdm.contrib.concurrent import process_map
from tqdm import tqdm
import numpy as np
from designer import MaxDense, AfReliedDesigner
import pandas as pd


def process_item(_args) -> tuple:
    # 独立处理每个item的函数
    _i, _df = _args
    _md = MaxDense(_df, 'E:\我的坚果云\ECS_mid')
    _md.form_chain(_i)
    return _i, _md.get_risk()


if __name__ == '__main__':
    os.chdir('E:\我的坚果云\ECS_final')
    df = pd.read_excel('variants.xlsx', index_col=0)
    path = 'E:\我的坚果云\ECS_mid'

    af = AfReliedDesigner(df, path)


    columns = ['n', 'md', 'af']
    df_result = pd.DataFrame(columns=columns)


    items = np.arange(1, 5001, 10)  # 待处理数据
    dfs = len(items) * [df]
    args = zip(items, dfs)
    results = process_map(process_item, args, max_workers=8, chunksize=10, smoothing=0.3)

    n_ = [results[i][0] for i in range(len(results))]
    mds = [results[i][1] for i in range(len(results))]
    afs = [af.result[n] for n in n_]

    df_result['n'] = n_
    df_result['md'] = mds
    df_result['af'] = afs

    df_result.to_excel('tmp.xlsx', index=True)

