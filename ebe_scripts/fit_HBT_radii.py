#!/usr/bin/env python3

import sys

import h5py
import numpy as np
from scipy.optimize import curve_fit

hbarc = 0.19733
HBT_PREFIX = "HBT_correlation_function_KT_"
RADII_PREFIX = "HBT_radii_KT_"
Q_CUTS = (0.05, 0.075, 0.1, 0.125, 0.15)
HEADER = (
    "# q_cut[GeV]  lambda  lambda_err  R_out[fm]  R_out_err[fm]  "
    "R_side[fm]  R_side_err[fm]  R_long[fm]  R_long_err[fm]  "
    "R_os[fm]  R_os_err[fm]  R_ol[fm]  R_ol_err[fm]"
)


def gaussian_3d(q_arr, lambda_, R_out, R_side, R_long, R_os, R_ol):
    """Evaluate the Gaussian parametrization from arXiv:1403.4972v1."""
    q_out, q_side, q_long = q_arr
    exponent = (
        (R_out * q_out)**2 + (R_side * q_side)**2 + (R_long * q_long)**2
        + 2 * q_out * q_side * R_os**2 + 2 * q_out * q_long * R_ol**2
    )
    return lambda_ * np.exp(-exponent / hbarc**2)


def fit_radii(data, label):
    if data.ndim != 2 or data.shape[1] < 8:
        raise ValueError(f"{label} requires the eight-column averaged HBT format")
    rows = len(data)
    nq = round(rows ** (1 / 3))
    if not rows or nq**3 != rows:
        raise ValueError(f"{label} does not have a cubic HBT grid ({rows} rows)")

    q = data[:, :3]
    correlation, error = data[:, 6:8].T
    valid = (np.all(np.isfinite(q), axis=1) & np.isfinite(correlation)
             & np.isfinite(error) & (error > 0))
    magnitude = np.linalg.norm(q, axis=1)
    output = np.full((len(Q_CUTS), 13), np.nan)
    output[:, 0] = Q_CUTS

    for index, q_cut in enumerate(Q_CUTS):
        selected = valid & (magnitude > 0) & (magnitude < q_cut)
        if np.count_nonzero(selected) <= 6:
            continue
        try:
            params, covariance = curve_fit(
                gaussian_3d, q[selected].T, correlation[selected],
                p0=(1.0, 5.0, 5.0, 5.0, 0.1, 0.1),
                sigma=error[selected], absolute_sigma=True,
            )
        except RuntimeError as error_message:
            print(f"{label}, q < {q_cut:g} GeV: {error_message}", file=sys.stderr)
            continue
        output[index, 1::2] = params
        output[index, 2::2] = np.sqrt(np.diag(covariance))
    return output


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} database_file.h5", file=sys.stderr)
        return 1
    with h5py.File(sys.argv[1], "r+") as database:
        for group in database.values():
            if not isinstance(group, h5py.Group):
                continue
            for key in list(group):
                if not key.startswith(HBT_PREFIX) or not key.endswith(".dat"):
                    continue
                label = f"{group.name}/{key}"
                print(f"Analyzing {label} ...")
                output = fit_radii(np.asarray(group[key]), label)
                output_name = RADII_PREFIX + key[len(HBT_PREFIX):]
                if output_name in group:
                    del group[output_name]
                result = group.create_dataset(
                    output_name, data=output, compression="gzip", compression_opts=9
                )
                result.attrs["header"] = HEADER
    return 0


if __name__ == "__main__":
    sys.exit(main())
