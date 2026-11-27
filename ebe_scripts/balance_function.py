import numpy as np


def normalized_correlation(same, mixed, same_error, mixed_error):
    """Return the normalized ratio and error, masking unsupported bins."""
    correlation = np.full(same.shape, np.nan)
    error = np.full(same.shape, np.nan)
    total_same, total_mixed = same.sum(), mixed.sum()
    if total_same > 0 and total_mixed > 0:
        valid = mixed > 0
        norm = total_mixed / total_same
        ratio = same[valid] / mixed[valid]
        correlation[valid] = norm * ratio
        error[valid] = norm * np.hypot(
            same_error[valid], ratio * mixed_error[valid]
        ) / mixed[valid]
    return correlation, error
