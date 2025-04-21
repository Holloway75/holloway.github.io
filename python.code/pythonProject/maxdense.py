import pandas as pd
import numpy as np
import os


def prob_merge(arr_:np.ndarray):
    return 1 - np.exp(np.sum(np.log(1 - arr_)))


class Preprocessor:
    def __init__(self, df:pd.DataFrame, ban_variants_list=None):
        if ban_variants_list:
            df.drop(ban_variants_list, inplace=True)
        df.sort_values(by='af', ascending=False, inplace=True)

        columns = ['af', 'cont', 'dense', 'gene', 'id']
        self.data = pd.DataFrame(columns=columns)
        self.data = self.data.astype({'af':'float64', 'cont':'float64', 'dense':'float64', 'gene':'object', 'id':'object'})

        self.data[['af', 'gene', 'id']] = df[['af', 'gene', 'id']]
        second_level = self.data.groupby('gene').cumcount()
        self.data.set_index(['gene', second_level], inplace=True)


class MaxDense(Preprocessor):
    def __init__(self, df, ban_variants_list=None):
        super().__init__(df=df, ban_variants_list=ban_variants_list)
        self.remain_loci = None
        self.panel_size = None
        self.result_genes = []
        self.result_varid = []
        self.result_counts = []
        self.result = None

        gene_list = list(set(df['gene']))
        record = {
            'starts': np.zeros(len(gene_list)),
            'var_counts': [self.data.loc[gene].shape[0] for gene in gene_list]
        }
        self.records = pd.DataFrame(data=record, index=gene_list)
        self.records = self.records.astype({'starts': 'int32', 'var_counts': 'int32'})
        self.get_cont()


    def get_cont(self):
        auto_list, xlink_list = self._get_inherit_type_list()
        for gene, var_counts in zip(self.records.index, self.records['var_counts']):
            afs = np.array(self.data.xs(gene)['af'])
            merged_afs = 2 * np.array([0] + [prob_merge(afs[:i+1]) for i in range(var_counts)])
            if gene in auto_list:
                merged_afs = merged_afs ** 2
            elif gene in xlink_list:
                pass
            else:
                raise ValueError()

            self.data.loc[gene, 'cont'] = [merged_afs[i+1] - merged_afs[i] for i in range(var_counts)]
            self.data.loc[gene, 'dense'] = self._get_dense(gene, 0, var_counts)


    @staticmethod
    def _get_inherit_type_list():
        os.chdir('E:\我的坚果云\ECS_mid')
        with open('gene.autosome.list', 'r') as f:
            auto_list = f.read().splitlines()
        with open('gene.x_link.list', 'r') as f:
            xlink_list = f.read().splitlines()
        return auto_list, xlink_list


    def _get_dense(self, gene, start, stop):
        return [np.mean(self.data.loc[gene, 'cont'][start:i+1]) for i in np.arange(start, stop)]


    def form_chain(self, panel_size):
        self.remain_loci = panel_size
        while self.remain_loci:
            self._tail_removal()
            self._prolong_chain()

        data = {
            'gene': self.result_genes,
            'varid': self.result_varid,
            'counts':self.result_counts
        }
        self.result = pd.DataFrame(data=data)

    def _tail_removal(self):
        condition = self.records['var_counts'] > self.remain_loci
        _indices = self.records[condition].index

        if len(_indices) == 0:
            pass
        else:
            for gene in _indices:
                _start = self.records.loc[gene, 'starts'] + self.remain_loci
                _end = self.records.loc[gene, 'starts'] + self.records.loc[gene, 'var_counts']
                col_data = self.data.loc[gene, 'dense'].values
                col_data[_start: _end] = 0
                self.data.loc[gene, 'dense'] = col_data
                self.records.loc[gene, 'var_counts'] = self.remain_loci


    def _prolong_chain(self):
        max_id = self.data['dense'].idxmax()
        gene = max_id[0]
        _start = self.records.loc[gene, 'starts']
        _var_counts = max_id[1] - _start + 1
        _end = _start + _var_counts

        self.records.loc[gene, 'starts'] += _var_counts
        self.result_genes.append(gene)
        self.result_varid.append(max_id[1])
        self.result_counts.append(_var_counts)
        self.remain_loci -= _var_counts

        col_data = self.data.loc[gene, 'dense'].values
        col_data[_start: _end] = 0
        self.data.loc[gene, 'dense'] = col_data

