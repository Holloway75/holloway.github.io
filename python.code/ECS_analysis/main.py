import os

import pandas as pd

import data_prepare
from addresss import *


if __name__ == '__main__':
   df = pd.read_csv('sample.combined.csv')

   df_area = data_prepare.convert_in_areas(df).set_index('area')
   print(data_prepare.transform_merge_area(df_area, Area_counterparts))
   print(pd.read_csv('areas.merged.csv', dtype='int8'))




