def multiplicity_bounds(descending_multiplicities, percentile_min, percentile_max):
    """Return (lower, upper) for selection with lower < multiplicity <= upper.

    Multiplicity ties at an internal edge stay together in the more peripheral
    bin. The last bin includes the lowest multiplicity.
    """
    count = len(descending_multiplicities)
    if not count or percentile_min >= percentile_max:
        return float("inf"), float("inf")

    upper = descending_multiplicities[int(count * percentile_min / 100)]
    lower = (
        float("-inf")
        if percentile_max == 100
        else descending_multiplicities[int(count * percentile_max / 100)]
    )
    return lower, upper
