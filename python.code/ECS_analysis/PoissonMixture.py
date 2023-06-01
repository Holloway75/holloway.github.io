import numpy as np
from sklearn import cluster
from scipy.special import logsumexp
import warnings

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
        nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
        self.weights = nk / self.n_samples
        self.params = np.dot(resp.T, X) / nk[:, np.newaxis]

    def _estimate_log_prob(self, X):
        par = self.params
        return -par.T + X * np.log(par.T) - np.array([np.log(np.arange(1, k+1)).sum() for k in X]).reshape(-1, 1)

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
