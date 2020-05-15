#!/usr/bin/env python3

from os import path
import sys
from numpy import *
from scipy.optimize import curve_fit

import h5py

hbarc = 0.19733
eps = 1e-15

def gaussian_3d(q_arr, lambda_, R_out, R_side, R_long, R_os, R_ol):
    """ the fit function is according to arXiv: 1403.4972v1 """
    (q_out, q_side, q_long) = q_arr

    R_out /= hbarc          # [1/GeV]
    R_side /= hbarc
    R_long /= hbarc
    R_os /= hbarc
    R_ol /= hbarc
    gauss = (lambda_*exp(
            - (  (R_out*q_out)**2 + (R_side*q_side)**2 + (R_long*q_long)**2
               + 2.*q_out*q_side*R_os**2. + 2.*q_out*q_long*R_ol**2.)
            )
    )

    # flatten the 3D Gaussian down to 1D
    return ravel(gauss)

if len(sys.argv) < 2:
    print("Usage: {} database_file.h5".format(sys.argv[0]))
    exit(0)

datafile_name = sys.argv[1]

q_cut_max_list = [0.05, 0.075, 0.1, 0.125, 0.15]
KT_cut_list = ['0_0.2', '0.2_0.4', '0.4_0.6', '0.6_0.8']

q_cut_min = 0.

header = (  "# q_cut[GeV]  lambda  lambda_err  R_out[fm]  R_out_err[fm]  "
          + "R_side[fm]  R_side_err[fm]  R_long [fm]  R_long_err[fm]  "
          + "R_os[fm]  R_os_err[fm]  R_ol[fm]  R_ol_err[fm]")

# load theory database
th_database = h5py.File(datafile_name)
event_list = list(th_database.keys())
for ievent in event_list:
    data_group = th_database[ievent]
    for KT_cut in KT_cut_list:
        filename="HBT_correlation_function_KT_{}.dat".format(KT_cut)

        HBTradii_filename = "HBT_radii_KT_{}.dat".format(KT_cut)
        if HBTradii_filename in data_group.keys():
            del data_group[HBTradii_filename]

        print("Analyzing {}: {} ...".format(ievent, filename))

        HBT_data = nan_to_num(data_group.get(filename))
        nq = int(HBT_data.shape[0]**(1/3))+1

        q_out = HBT_data[:, 0].reshape(nq, nq, nq)
        q_side = HBT_data[:, 1].reshape(nq, nq, nq)
        q_long = HBT_data[:, 2].reshape(nq, nq, nq)
        HBT_Corr = HBT_data[:, 6].reshape(nq, nq, nq)
        HBT_Corr_err = HBT_data[:, 7].reshape(nq, nq, nq) + eps

        output = []
        for q_cut in q_cut_max_list:
            idx = (   (sqrt(q_out**2. + q_side**2. + q_long**2.) > q_cut_min)
                    & (sqrt(q_out**2. + q_side**2. + q_long**2.) < q_cut))

            q_arr = [q_out[idx], q_side[idx], q_long[idx]]
            guess_vals = [1.0, 5., 5., 5., 0.1, 0.1]

            fit_params, cov_mat = curve_fit(
                    gaussian_3d, q_arr, ravel(HBT_Corr[idx]), p0=guess_vals,
                    sigma=ravel(HBT_Corr_err[idx]), absolute_sigma=True)
            fit_errors = sqrt(diag(cov_mat))
            fit_corr = gaussian_3d(q_arr, *fit_params).reshape(
                                                        HBT_Corr[idx].shape)

            # manually calculate R-squared goodness of fit
            fit_residual = HBT_Corr[idx] - fit_corr
            fit_Rsquared = 1 - var(fit_residual)/var(HBT_Corr[idx])

            temp = []
            for x, y in zip(fit_params, fit_errors):
                temp += [x, y]

            output.append([q_cut] + temp)

            if q_cut == 0.1:
                print("R_long = {} fm, R_side = {} fm, R_long = {} fm".format(
                    fit_params[1], fit_params[2], fit_params[3]))

        h5data = data_group.create_dataset(
                    "{0}".format(HBTradii_filename), data=array(output),
                    compression="gzip", compression_opts=9)
        h5data.attrs.create("header", string_(header))
th_database.close()

