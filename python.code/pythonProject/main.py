import os

import numpy as np
import pandas as pd

if __name__ == '__main__':
    os.chdir('E:\我的坚果云\ECS_final')
    df = pd.read_excel('variants.xlsx', index_col=0)
    a1 = np.array((1,2,3))
    a2 = np.array((2,3,4))
    print(a1 * a2)