from unittest.mock import inplace

import pandas as pd
import numpy as np
import os

from fontTools.merge.util import first


def prob_merge(arr_:np.ndarray):
    return 1 - np.exp(np.sum(np.log(1 - arr_)))


class Preprocessor:
    def __init__(self, df:pd.DataFrame, ban_variants_list=None):
        """

        :param df: 'af', 'gene', 'inheritance' were required in columns.
        :param ban_variants_list:
        """
        if ban_variants_list:
            df.drop(ban_variants_list, inplace=True)
        df.sort_values(by='af', ascending=False, inplace=True)
        self.gene_list = list(set(df['gene']))

        columns = ['af', 'cont', 'dense', 'gene', 'id']
        self.data = pd.DataFrame(columns=columns)
        self.data = self.data.astype({'af':'float64', 'cont':'float64', 'dense':'float64', 'gene':'object', 'id':'object'})

        self.data[['af', 'gene', 'id']] = df[['af', 'gene', 'id']]
        second_level = self.data.groupby('gene').cumcount()
        self.data.set_index(['gene', second_level], inplace=True)


class MaxDense(Preprocessor):
    def __init__(self, df, ban_variants_list=None):
        super().__init__(df=df, ban_variants_list=ban_variants_list)
        self.get_cont()

    def get_cont(self):
        auto_list, xlink_list = self._get_inherit_type_list()
        for gene in self.gene_list:



            afs = np.array(self.data.xs(gene)['af'])
            var_counts = len(afs)
            merged_afs = np.array([0] + [prob_merge(afs[:i+1]) for i in range(var_counts)])
            if gene in auto_list:
                merged_afs = merged_afs ** 2
            elif gene in xlink_list:
                pass
            else:
                raise ValueError()

            self.data.loc[gene, 'cont'] = [merged_afs[i+1] - merged_afs[i] for i in range(var_counts)]

            print(self.data.loc[gene].dtypes)
            print(self.data.loc[gene])
            exit()


    @staticmethod
    def _get_inherit_type_list():
        os.chdir('E:\我的坚果云\ECS_mid')
        with open('gene.autosome.list', 'r') as f:
            auto_list = f.read().splitlines()
        with open('gene.x_link.list', 'r') as f:
            xlink_list = f.read().splitlines()
        return auto_list, xlink_list