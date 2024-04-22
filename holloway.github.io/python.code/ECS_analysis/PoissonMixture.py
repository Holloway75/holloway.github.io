import numpy as np
import pandas as pd
from sklearn import cluster
from scipy.special import logsumexp
from addresss import *
import warnings
from scipy.spatial import KDTree

__all__ = [
    'PoissonMixture',
    'IndividualRiskEstimator'
]


class PoissonMixture:
    """
    通过EM算法计算泊松混合模型的参数
    """
    def __init__(
            self,
            n_components=1,
            *,
            tol=1e-3,
            max_iter=200,
    ):
        self.n_components = n_components
        self.tol = tol
        self.max_iter = max_iter
        self.converged_ = False

    def _initialize_parameters(self, X):
        self.n_samples = X.shape[0]
        resp = np.zeros((self.n_samples, self.n_components))
        label = (
            cluster.KMeans(
                n_clusters=self.n_components, n_init=1
            )
            .fit(X)
            .labels_
        )
        resp[np.arange(self.n_samples), label] = 1
        nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
        self.weights = nk / self.n_samples
        self.params = np.dot(resp.T, X) / nk[:, np.newaxis]

    def fit(self, X):
        self._initialize_parameters(X)
        lower_bound = -np.inf
        for n_iter in range(1, self.max_iter + 1):
            prev_lower_bound = lower_bound

            log_prob_norm, log_resp = self._e_step(X)
            self._m_step(X, log_resp)

            lower_bound = log_prob_norm
            change = lower_bound - prev_lower_bound
            if abs(change) < self.tol:
                self.converged_ = True
                break

        if not self.converged_ and self.max_iter > 0:
            warnings.warn(
                "Initialization did not converge. "
                "Try different init parameters, "
                "or increase max_iter, tol "
                "or check for degenerate data.",
                DeprecationWarning,
            )
        return self

    def _e_step(self, X):
        log_prob_norm, log_resp = self._estimate_log_prob_resp(X)
        return np.mean(log_prob_norm), log_resp

    def _m_step(self, X, log_resp):
        resp = np.exp(log_resp)
        nk = resp.sum(axis=0) + 10*np.finfo(resp.dtype).eps
        self.weights = nk / self.n_samples
        self.params = np.dot(resp.T, X) / nk[:, np.newaxis]

    def _estimate_log_prob(self, X):
        par = self.params
        eps = np.finfo(par.dtype).eps
        return -par.T + X * np.log(par.T+eps) - np.array([np.log(np.arange(1, k+1)).sum() for k in X]).reshape(-1, 1)

    def _estimate_log_weights(self):
        return np.log(self.weights)

    def _estimate_weighted_log_prob(self, X):
        return self._estimate_log_prob(X) + self._estimate_log_weights()

    def _estimate_log_prob_resp(self, X):
        weighted_log_prob = self._estimate_weighted_log_prob(X)
        log_prob_norm = logsumexp(weighted_log_prob, axis=1)
        with np.errstate(under="ignore"):
            # ignore underflow
            log_resp = weighted_log_prob - log_prob_norm[:, np.newaxis]
        return log_prob_norm, log_resp

    def score(self, X):
        log_prob_sum = logsumexp(self._estimate_weighted_log_prob(X), axis=1)
        return log_prob_sum.mean()

    def bic(self, X):
        k = self.n_components + self.n_components - 1
        return -2 * self.score(X) * X.shape[0] + k * np.log(X.shape[0])


class IndividualRiskEstimator:
    """
    估计指定基因的个体化风险
    """
    def __init__(self, input_df: pd.DataFrame, xlink=False, pc=1):
        if xlink == True:
            self.data = input_df[input_df['sex']==0].sort_values(by='PC1').reset_index(drop=False)
            self.glist = Xlink_list
        else:
            self.data = input_df.sort_values(by='PC1').reset_index(drop=False)
            self.glist = Auto_list
        self._sample_counts = self.data.shape[0]

    def estimate_indiv_risk(self, pc):
        return self._get_indiv_risk(self.glist, pc)

    def carrier_matrix(self):
        car_mat = self._carrier_matrix(self.glist)
        return pd.DataFrame(columns=self.glist, index=self.data['ecs_id'], data=car_mat)

    def _carrier_matrix(self, gene_list):
        car_mat = np.zeros((self._sample_counts, len(gene_list)), dtype='i1')
        df_tmp = self.data.dropna(axis=0, subset='gene')
        for gene in gene_list:
            labels = [i for i in df_tmp.index if gene in df_tmp.loc[i, 'gene'].split(":")]
            car_mat[labels, [gene_list.index(gene)]*len(labels)] = 1
        return car_mat

    def _get_neighbors(self, id, pc):
        column = ['PC%d' % i for i in np.arange(1, pc+1)]
        rng = KDTree(data=np.array(self.data[column]).reshape(-1, pc))
        _, neighbor_index = rng.query(self.data.loc[id, column], 1000)
        return neighbor_index

    def _get_indiv_risk(self, gene_list, pc):
        car_mat = self._carrier_matrix(gene_list)
        indiv_mat = np.zeros((self._sample_counts - 1000, car_mat.shape[1]), dtype='i2')
        index_list = self.data.index.tolist()[500:-500]
        t = 0
        for i in index_list:
            indiv_mat[t] = np.sum(car_mat[self._get_neighbors(i, pc)], axis=0)
            t += 1
        return pd.DataFrame(data=indiv_mat, columns=gene_list, index=self.data['ecs_id'].tolist()[500:-500])




