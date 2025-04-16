import datetime

import pandas as pd
import numpy as np
from statsmodels.sandbox.distributions.genpareto import shape
from statsmodels.sandbox.stats.multicomp import TukeyHSDResults


def prob_merge(arr_:np.ndarray):
    return 1 - np.exp(np.sum(np.log(1 - arr_)))


class Preprocessor:
    def __init__(self, df:pd.DataFrame, ban_variants_list=None):
        """

        :param df: 'af', 'gene', 'inheritance' were required in columns.
        :param ban_variants_list:
        """
        df.sort_values(by='gene', inplace=True)
        if ban_variants_list:
            df.drop(ban_variants_list, inplace=True)

        columns = ['af', 'cont', 'dense', 'gene']
        self.data = pd.DataFrame(columns=columns, index=df.index)
        self.data[['af', 'gene']] = df[['af', 'gene']]


    def get_cont(self, dict_gene_risk):
        pass






class MaxDense(Preprocessor):
    def __init__(self, df, ban_variants_list=None):
        super().__init__(df=df, ban_variants_list=ban_variants_list)
        self.dict_cont_perct = dict()


    def _get_cont(self, gene):
        pass

    def get_dense(self):
        pass