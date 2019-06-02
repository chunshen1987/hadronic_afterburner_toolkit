#! /usr/bin/env python
"""
     This script performs event averaging for particle 
     spectra and anisotropic flow coefficients calculated 
     from event-by-event simulations

     v_n is analyzed up to n = 6

     Format for particle_XXX_vndata.dat file:
     n_order  real_part  real_part_err  imag_part  imag_part_err
     
     Format for particle_XXX_vndata_diff.dat file:
     pT(GeV)  pT_err(GeV)  dN/(2pi dy pT dpT)(GeV^-2)  dN/(2pi dy pT dpT)_err(GeV^-2)
     vn_real  vn_real_err  vn_imag  vn_imag_err

     All the errors are only statistic errors
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
from numpy import *
import h5py
import shutil

try:
    data_path = path.abspath(argv[1])
    data_name = data_path.split("/")[-1]
    results_folder_name = data_name.split(".h5")[0]
    avg_folder = path.join(path.abspath(argv[2]),
                           results_folder_name)
    print("output folder: %s" % avg_folder)
    if(path.isdir(avg_folder)):
        print("folder %s already exists!" % avg_folder)
        var = input("do you want to delete it? [y/N]")
        if 'y' in var:
            shutil.rmtree(avg_folder)
        else:
            print("please choose another folder path~")
            exit(0)
    mkdir(avg_folder)
except IndexError:
    print("Usage: {} results.h5 results_folder".format(str(argv[0])))
    exit(1)

ecc_filename = "eccentricities_evo_eta_-0.5_0.5.dat"
ep_filename  = "momentum_anisotropy_eta_-0.5_0.5.dat"

tau = linspace(0.4, 5.4, 1001)
ntau = len(tau)
hf = h5py.File(data_path, "r")
event_list = list(hf.keys())
nev = len(event_list)

ep2 = zeros([ntau, 7])
ep3 = zeros([ntau, 7])
ecc2_dot_ep2 = zeros([ntau, 7])
ecc3_dot_ep3 = zeros([ntau, 7])
ep2[:, 0] = tau
ep3[:, 0] = tau
ecc2_dot_ep2[:, 0] = tau
ecc3_dot_ep3[:, 0] = tau

for ifolder, event_name in enumerate(event_list):
    event_group = hf.get(event_name)
    ecc_n_data = array(event_group.get(ecc_filename))
    ep_n_data  = array(event_group.get(ep_filename))

    ntau_loc = len(ecc_n_data[:, 0])

    ep_2_ideal = sqrt(ep_n_data[:, 4]**2. + ep_n_data[:, 5]**2.)
    ep_2_shear = sqrt(ep_n_data[:, 6]**2. + ep_n_data[:, 7]**2.)
    ep_2_full = sqrt(ep_n_data[:, 8]**2. + ep_n_data[:, 9]**2.)
    ep_3_ideal = sqrt(ep_n_data[:, 10]**2. + ep_n_data[:, 11]**2.)
    ep_3_shear = sqrt(ep_n_data[:, 12]**2. + ep_n_data[:, 13]**2.)
    ep_3_full = sqrt(ep_n_data[:, 14]**2. + ep_n_data[:, 15]**2.)
    ecc2_dot_ep2_ideal = (
        ecc_n_data[:, 3]*ep_n_data[:, 4] + ecc_n_data[:, 4]*ep_n_data[:, 5])
    ecc2_dot_ep2_shear = (
        ecc_n_data[:, 3]*ep_n_data[:, 6] + ecc_n_data[:, 4]*ep_n_data[:, 7])
    ecc2_dot_ep2_full = (
        ecc_n_data[:, 3]*ep_n_data[:, 8] + ecc_n_data[:, 4]*ep_n_data[:, 9])
    ecc3_dot_ep3_ideal = (
        ecc_n_data[:, 5]*ep_n_data[:, 10] + ecc_n_data[:, 6]*ep_n_data[:, 11])
    ecc3_dot_ep3_shear = (
        ecc_n_data[:, 5]*ep_n_data[:, 12] + ecc_n_data[:, 6]*ep_n_data[:, 13])
    ecc3_dot_ep3_full = (
        ecc_n_data[:, 5]*ep_n_data[:, 14] + ecc_n_data[:, 6]*ep_n_data[:, 15])

    if ntau_loc < ntau:
        ep2[:ntau_loc, 1] += ep_2_ideal
        ep2[:ntau_loc, 2] += ep_2_ideal**2.
        ep2[:ntau_loc, 3] += ep_2_shear
        ep2[:ntau_loc, 4] += ep_2_shear**2.
        ep2[:ntau_loc, 5] += ep_2_full
        ep2[:ntau_loc, 6] += ep_2_full**2.
        ep3[:ntau_loc, 1] += ep_3_ideal
        ep3[:ntau_loc, 2] += ep_3_ideal**2.
        ep3[:ntau_loc, 3] += ep_3_shear
        ep3[:ntau_loc, 4] += ep_3_shear**2.
        ep3[:ntau_loc, 5] += ep_3_full
        ep3[:ntau_loc, 6] += ep_3_full**2.
        ecc2_dot_ep2[:ntau_loc, 1] += ecc2_dot_ep2_ideal
        ecc2_dot_ep2[:ntau_loc, 2] += ecc2_dot_ep2_ideal**2.
        ecc2_dot_ep2[:ntau_loc, 3] += ecc2_dot_ep2_shear
        ecc2_dot_ep2[:ntau_loc, 4] += ecc2_dot_ep2_shear**2.
        ecc2_dot_ep2[:ntau_loc, 5] += ecc2_dot_ep2_full
        ecc2_dot_ep2[:ntau_loc, 6] += ecc2_dot_ep2_full**2.
        ecc3_dot_ep3[:ntau_loc, 1] += ecc3_dot_ep3_ideal
        ecc3_dot_ep3[:ntau_loc, 2] += ecc3_dot_ep3_ideal**2.
        ecc3_dot_ep3[:ntau_loc, 3] += ecc3_dot_ep3_shear
        ecc3_dot_ep3[:ntau_loc, 4] += ecc3_dot_ep3_shear**2.
        ecc3_dot_ep3[:ntau_loc, 5] += ecc3_dot_ep3_full
        ecc3_dot_ep3[:ntau_loc, 6] += ecc3_dot_ep3_full**2.
    else:
        ep2[:, 1] += ep_2_ideal[:ntau]
        ep2[:, 2] += ep_2_ideal[:ntau]**2.
        ep2[:, 3] += ep_2_shear[:ntau]
        ep2[:, 4] += ep_2_shear[:ntau]**2.
        ep2[:, 5] += ep_2_full[:ntau]
        ep2[:, 6] += ep_2_full[:ntau]**2.
        ep3[:, 1] += ep_3_ideal[:ntau]
        ep3[:, 2] += ep_3_ideal[:ntau]**2.
        ep3[:, 3] += ep_3_shear[:ntau]
        ep3[:, 4] += ep_3_shear[:ntau]**2.
        ep3[:, 5] += ep_3_full[:ntau]
        ep3[:, 6] += ep_3_full[:ntau]**2.
        ecc2_dot_ep2[:, 1] += ecc2_dot_ep2_ideal[:ntau]
        ecc2_dot_ep2[:, 2] += ecc2_dot_ep2_ideal[:ntau]**2.
        ecc2_dot_ep2[:, 3] += ecc2_dot_ep2_shear[:ntau]
        ecc2_dot_ep2[:, 4] += ecc2_dot_ep2_shear[:ntau]**2.
        ecc2_dot_ep2[:, 5] += ecc2_dot_ep2_full[:ntau]
        ecc2_dot_ep2[:, 6] += ecc2_dot_ep2_full[:ntau]**2.
        ecc3_dot_ep3[:, 1] += ecc3_dot_ep3_ideal[:ntau]
        ecc3_dot_ep3[:, 2] += ecc3_dot_ep3_ideal[:ntau]**2.
        ecc3_dot_ep3[:, 3] += ecc3_dot_ep3_shear[:ntau]
        ecc3_dot_ep3[:, 4] += ecc3_dot_ep3_shear[:ntau]**2.
        ecc3_dot_ep3[:, 5] += ecc3_dot_ep3_full[:ntau]
        ecc3_dot_ep3[:, 6] += ecc3_dot_ep3_full[:ntau]**2.

ep2[:, 1:] /= nev
ep2[:, 2::2] = sqrt(ep2[:, 2::2] - ep2[:, 1::2]**2.)/sqrt(nev)
ep3[:, 1:] /= nev
ep3[:, 2::2] = sqrt(ep3[:, 2::2] - ep3[:, 1::2]**2.)/sqrt(nev)
ecc2_dot_ep2[:, 1:] /= nev
ecc2_dot_ep2[:, 2::2] = (
        sqrt(ecc2_dot_ep2[:, 2::2] - ecc2_dot_ep2[:, 1::2]**2.)/sqrt(nev))
ecc3_dot_ep3[:, 1:] /= nev
ecc3_dot_ep3[:, 2::2] = (
        sqrt(ecc3_dot_ep3[:, 2::2] - ecc3_dot_ep3[:, 1::2]**2.)/sqrt(nev))

filename = 'momentum_aniso_ep2_evo.dat'
header_text = ("# tau  ep2_ideal  ep2_ideal_err  ep2_shear  ep2_shear_err  "
                + "ep2_full  ep2_full_err")
savetxt(filename, ep2, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'momentum_aniso_ep3_evo.dat'
header_text = ("# tau  ep3_ideal  ep3_ideal_err  ep3_shear  ep3_shear_err  "
                + "ep3_full  ep3_full_err")
savetxt(filename, ep3, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'ecc2_dot_ep2_evo.dat'
header_text = ("# tau  ideal  ideal_err  shear  shear_err  full  full_err")
savetxt(filename, ecc2_dot_ep2, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'ecc3_dot_ep3_evo.dat'
header_text = ("# tau  ideal  ideal_err  shear  shear_err  full  full_err")
savetxt(filename, ecc3_dot_ep3, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)
