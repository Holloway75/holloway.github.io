import os
from maxdense import MaxDense
import numpy as np
import pandas as pd

if __name__ == '__main__':
    os.chdir('E:\我的坚果云\ECS_final')
    df = pd.read_excel('variants.xlsx', index_col=0)
    md = MaxDense(df)



