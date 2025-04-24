import pandas as pd
import numpy as np
import os


def prob_merge(arr_:np.ndarray):
    """Calculate merged probability using complementary log-log model.

    Args:
        arr_: Array of probabilities to merge

    Returns:
        Merged probability calculated as 1 - product(1 - p) for p in arr_
    """
    return 1 - np.exp(np.sum(np.log(1 - arr_)))


class BaseDesigner:
    """Base class for panel design optimization strategies with path flexibility."""

    def __init__(self, df: pd.DataFrame,
                 gene_list_path: str,
                 ban_variants_list=None):
        """Initialize base designer with configurable gene list path.

        Args:
            df: Input DataFrame containing variant data
            gene_list_path: Customizable path to inheritance type lists
            ban_variants_list: Variants to exclude from analysis
        """
        self.gene_list_path = gene_list_path
        if ban_variants_list:
            df.drop(ban_variants_list, inplace=True)
        df.sort_values(by='af', ascending=False, inplace=True)
        self.data = None
        self.records = None
        self.result = None
        self.merged_risks = dict()   # Stores merged risk values per gene
        self._make_data(df)
        self._merge_risks()


    def _make_data(self, df) -> None:
        """Organize data into hierarchical structure (gene, variant order)."""
        columns = ['af', 'gene', 'id']
        self.data = pd.DataFrame(columns=columns)
        self.data = self.data.astype({'af': 'float64', 'gene': 'object', 'id': 'object'})

        # Create multi-index: (gene, variant order within gene)
        self.data[['af', 'gene', 'id']] = df[['af', 'gene', 'id']]
        second_level = self.data.groupby('gene').cumcount()
        self.data.set_index(['gene', second_level], inplace=True)

        # Initialize gene tracking records
        gene_list = list(set(df['gene']))
        record = {
            'starts': np.zeros(len(gene_list)),     # Current selection pointer
            'var_counts': [self.data.loc[gene].shape[0] for gene in gene_list]  # Total variants
        }
        self.records = pd.DataFrame(data=record, index=gene_list)
        self.records = self.records.astype({'starts': 'int32', 'var_counts': 'int32'})


    def _merge_risks(self) -> None:
        """Calculate cumulative risks using vectorized operations and inheritance types."""
        auto_list, xlink_list = self._get_inherit_type_list()
        for gene, var_counts in zip(self.records.index, self.records['var_counts']):
            gene_data = self.data.loc[gene]
            afs = gene_data['af'].values

            # Vectorized calculation of merged AF
            cumulative_log = np.log(1 - afs).cumsum()
            merged_afs = 2 * (1 - np.exp(cumulative_log))
            merged_afs = np.insert(merged_afs, 0, 0)

            # Apply inheritance pattern
            if gene in auto_list:
                merged_afs = merged_afs ** 2    # Autosomal recessive
            elif gene in xlink_list:
                pass    # X-linked (no squaring)
            else:
                raise ValueError(
                    f"Gene {gene} missing inheritance type annotation. "
                    f"Valid types: {auto_list + xlink_list}"
                )
            self.merged_risks[gene] = merged_afs


    def _get_inherit_type_list(self) -> tuple:
        """Load inheritance type lists from configurable path.

        Returns:
            Tuple containing (autosomal_genes, xlinked_genes)
        """
        auto_path = os.path.join(self.gene_list_path, 'gene.autosome.list')
        xlink_path = os.path.join(self.gene_list_path, 'gene.x_link.list')
        with open(auto_path, 'r') as f:
            auto_list = f.read().splitlines()
        with open(xlink_path, 'r') as f:
            xlink_list = f.read().splitlines()
        return auto_list, xlink_list


    def get_risk(self) -> float:
        """Calculate current total risk probability using merged values."""
        risks = [self.merged_risks[gene][_id] for gene, _id in zip(self.records.index, self.records['starts'])]
        return prob_merge(np.array(risks))


class MaxDense(BaseDesigner):
    """Implements optimized maximum density algorithm with incremental calculation."""

    def __init__(self, df: pd.DataFrame,
                 gene_list_path: str,
                 ban_variants_list=None):
        super().__init__(df=df, ban_variants_list=ban_variants_list, gene_list_path=gene_list_path)
        self.data['cont'] = pd.Series(dtype='float64')  # Marginal contribution
        self.data['dense'] = pd.Series(dtype='float64') # Average contribution
        self.remain_loci = None     # Remaining panel capacity
        self.panel_size = None      # Target panel size
        self.get_cont_dense()

        # Lightweight state preservation
        self.initial_records = self.records.to_dict()   # Stored as dictionary
        self.initial_dense = self.data['dense'].values.copy()    # Stored as numpy array


    def get_cont_dense(self) -> None:
        """Calculate contributions and densities using GeneChain objects."""
        for gene, var_counts in zip(self.records.index, self.records['var_counts']):
            # Calculate marginal contributions
            self.data.loc[gene, 'cont'] = [self.merged_risks[gene][i + 1] - self.merged_risks[gene][i] for i in
                                           range(var_counts)]

            # Incremental density calculation
            chain = GeneChain(self.data.loc[gene, 'cont'].tolist())
            self.data.loc[gene, 'dense'] = [chain.add_variant() for i in chain.cont]


    def form_chain(self, panel_size) -> None:
        """Construct optimal panel with length constraints.

        Args:
            panel_size: Target number of variants to include
        """
        self.remain_loci = panel_size
        result_genes, result_ids, result_counts, = [], [], []

        while self.remain_loci:
            self._tail_removal()        # Enforce length constraints
            self._prolong_chain(result_genes, result_ids, result_counts)    # Add optimal segment

        # Compile final panel
        data = {
            'gene': result_genes,
            'var_id': result_ids,
            'counts': result_counts,
        }
        self.result = pd.DataFrame(data=data)


    def _tail_removal(self) -> None:
        """Prune gene chains exceeding remaining capacity."""
        condition = self.records['var_counts'] > self.remain_loci
        _indices = self.records[condition].index

        if not _indices.empty:
            for gene in _indices:
                # Calculate non-viable segment range
                _start = self.records.loc[gene, 'starts'] + self.remain_loci
                _end = self.records.loc[gene, 'starts'] + self.records.loc[gene, 'var_counts']

                # Nullify non-viable densities
                col_data = self.data.loc[gene, 'dense'].values
                col_data[_start: _end] = 0
                self.data.loc[gene, 'dense'] = col_data
                self.records.loc[gene, 'var_counts'] = self.remain_loci


    def _prolong_chain(self,result_genes, result_ids, result_counts) -> None:
        """Select and integrate highest density segment."""
        max_id = self.data['dense'].idxmax()
        gene = max_id[0]
        _start = self.records.loc[gene, 'starts']
        _var_counts = max_id[1] - _start + 1
        _end = _start + _var_counts

        # Update tracking state
        self.records.loc[gene, 'starts'] += _var_counts
        self.records.loc[gene, 'var_counts'] -= _var_counts
        result_genes.append(gene)
        result_ids.append(max_id[1])
        result_counts.append(_var_counts)
        self.remain_loci -= _var_counts

        # Invalidate used segment
        col_data = self.data.loc[gene, 'dense'].values
        col_data[_start: _end] = 0
        self.data.loc[gene, 'dense'] = col_data


    def re_set(self) -> None:
        """Reset to initial state using lightweight preservation."""
        self.records = pd.DataFrame(self.initial_records)
        self.data['dense'] = self.initial_dense.copy()


class AfReliedDesigner(BaseDesigner):
    """Implements AF-based algorithm with risk progression tracking."""

    def __init__(self, df: pd.DataFrame,
                 gene_list_path: str,
                 ban_variants_list=None):
        super().__init__(df=df, gene_list_path=gene_list_path, ban_variants_list=ban_variants_list)
        self.data.set_index('af', append=True, inplace=True)
        self.data.sort_index(level=[2, 1], ascending=[False, True], inplace=True)
        self.result = dict()
        self._get_risk()


    def _get_risk(self) -> None:
        """Track risk changes when adding variants in AF order."""
        for num, index in enumerate(self.data.index):
            self.records.loc[index[0], 'starts'] = index[1] + 1
            self.result[num + 1] = self.get_risk()


class GeneChain:
    """Helper class for incremental density calculation."""

    def __init__(self, cont):
        """
        Args:
            cont: List of contribution values for a gene
        """
        self.cont = cont
        self.current_sum = 0.0
        self.current_length = 0

    def add_variant(self) -> float:
        """Incrementally add variant and return current density.

        Returns:
            Current average contribution (sum / length)
        """
        self.current_sum += self.cont[self.current_length]
        self.current_length += 1
        return self.current_sum / self.current_length

