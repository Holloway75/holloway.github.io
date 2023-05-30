import copy
from scipy.special import factorial
from scipy.stats import poisson
import numpy as np
from sklearn import cluster


__all__ = [
    'Poisson_Mixture'
]


class Poisson_Mixture:
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
        self.resp = resp
        nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
        self.weights = nk / self.n_samples
        self.params = np.dot(resp.T, X) / nk[:, np.newaxis]

    def fit(self, X):
        self._initialize_parameters(X)
        for iter in range(self.max_iter):
            self._e_step(X)
            self._m_step(X)
            print(self.params)
            if np.linalg.norm((self.params - self._pre_params)) < self.tol:
                break
        return self

    def _e_step(self, X):
        for i in np.arange(self.n_samples):
            pos_list = [poisson.pmf(X[i, 0], self.params[j, 0]) for j in range(self.n_components)]
            self.resp[i] = [j/sum(pos_list) for j in pos_list]

    def _m_step(self, X):
        self._pre_params = copy.deepcopy(self.params)
        nk = self.resp.sum(axis=0) + 10 * np.finfo(self.resp.dtype).eps
        self.weights = nk / self.n_samples
        self.params = np.dot(self.resp.T, X) / nk[:, np.newaxis]

    def bic(self, X):
        k = self.n_components + self.n_components - 1
        return -2 * self.score(X) * X.shape[0] + k * np.log(X.shape[0])


