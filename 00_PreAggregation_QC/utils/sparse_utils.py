import numpy as np
from scipy.sparse import issparse


def sparse_sum_rows(X) -> np.ndarray:
    """Sum all rows of a sparse or dense matrix. Returns 1-D float64 array of length n_cols."""
    if issparse(X):
        return np.asarray(X.sum(axis=0)).ravel().astype(np.float64)
    return np.asarray(X).sum(axis=0).astype(np.float64)


def donor_detection_vector(X) -> np.ndarray:
    """
    Return a boolean 1-D array of length n_cols indicating which columns
    have at least one non-zero entry across all rows of X.
    Used for prevalence counting (one donor at a time).
    """
    if issparse(X):
        return np.asarray((X > 0).sum(axis=0) > 0).ravel()
    return (np.asarray(X) > 0).any(axis=0)
