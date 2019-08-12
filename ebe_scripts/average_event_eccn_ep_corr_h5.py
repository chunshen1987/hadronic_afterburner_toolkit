#! /usr/bin/env python
"""
     This script computes the correlation between spatial eccentriicty
     and momentum anisotropy
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
from numpy import *
import h5py
import matplotlib.pyplot as plt
import shutil

def compute_psi_n(ep_n_cos, ep_n_sin, ecc_n_cos_init, ecc_n_sin_init):
    """Compute the angle between ep_n and ecc_n"""
    eccn_init_dot_epn_cos = ecc_n_cos_init*ep_n_cos + ecc_n_sin_init*ep_n_sin
    eccn_init_dot_epn_sin = ecc_n_cos_init*ep_n_sin - ecc_n_sin_init*ep_n_cos
    psi_n = arctan2(eccn_init_dot_epn_sin, eccn_init_dot_epn_cos)
    return psi_n

def pack_arrays(data, data_ideal, data_shear, data_full, ntau, ntau_loc):
    """package the first ntau elements to data"""
    if ntau_loc < ntau:
        data[:ntau_loc, 1] += data_ideal**2.
        data[:ntau_loc, 2] += data_ideal**4.
        data[:ntau_loc, 3] += data_shear**2.
        data[:ntau_loc, 4] += data_shear**4.
        data[:ntau_loc, 5] += data_full**2.
        data[:ntau_loc, 6] += data_full**4.
    else:
        data[:, 1] += data_ideal[:ntau]**2.
        data[:, 2] += data_ideal[:ntau]**4.
        data[:, 3] += data_shear[:ntau]**2.
        data[:, 4] += data_shear[:ntau]**4.
        data[:, 5] += data_full[:ntau]**2.
        data[:, 6] += data_full[:ntau]**4.

def compute_mean_and_stat_error(data, nev):
    """compute the mean and statistical error for data"""
    data[:, 1:] /= nev
    data[:, 2::2] = sqrt(data[:, 2::2] - data[:, 1::2]**2.)/sqrt(nev)
    data[:, 1::2] = sqrt(data[:, 1::2])
    data[:, 2::2] = data[:, 2::2]/(2.*data[:, 1::2] + 1e-16)

try:
    data_path = path.abspath(argv[1])
    data_name = data_path.split("/")[-1]
    results_folder_name = data_name.split(".h5")[0]
    avg_folder = path.join(path.abspath(argv[2]),
                           results_folder_name)
    print("output folder: %s" % avg_folder)
    if path.isdir(avg_folder):
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
ep_filename = "momentum_anisotropy_eta_-0.5_0.5.dat"

tau = linspace(0.4, 5.4, 1001)
#tau = linspace(0.4, 8.4, 1601)
ntau = len(tau)
hf = h5py.File(data_path, "r")
event_list = list(hf.keys())
nev = len(event_list)

ep2 = zeros([ntau, 7])
ep2_org = zeros([ntau, 7])
ep3 = zeros([ntau, 7])
psi_2 = zeros([ntau, 7])
psi_2_org = zeros([ntau, 7])
psi_3 = zeros([ntau, 7])
ecc2_dot_ep2 = zeros([ntau, 7])
ecc2_dot_ep2_org = zeros([ntau, 7])
ecc3_dot_ep3 = zeros([ntau, 7])
ep2[:, 0] = tau
ep2_org[:, 0] = tau
ep3[:, 0] = tau
psi_2[:, 0] = tau
psi_2_org[:, 0] = tau
psi_3[:, 0] = tau
ecc2_dot_ep2[:, 0] = tau
ecc2_dot_ep2_org[:, 0] = tau
ecc3_dot_ep3[:, 0] = tau

psi_2_init_list = []
psi_2_org_init_list = []
psi_3_init_list = []
psi_2_final_list = []
psi_2_org_final_list = []
psi_3_final_list = []

for ifolder, event_name in enumerate(event_list):
    event_group = hf.get(event_name)
    ecc_n_data = array(event_group.get(ecc_filename))
    ep_n_data = array(event_group.get(ep_filename))

    ntau_loc = len(ecc_n_data[:, 0])

    ep_2_org_ideal = sqrt(ep_n_data[:, 1]**2. + ep_n_data[:, 2]**2.)
    ep_2_org_shear = sqrt(ep_n_data[:, 3]**2. + ep_n_data[:, 4]**2.)
    ep_2_org_full = sqrt(ep_n_data[:, 5]**2. + ep_n_data[:, 6]**2.)
    ep_2_ideal = sqrt(ep_n_data[:, 7]**2. + ep_n_data[:, 8]**2.)
    ep_2_shear = sqrt(ep_n_data[:, 9]**2. + ep_n_data[:, 10]**2.)
    ep_2_full = sqrt(ep_n_data[:, 11]**2. + ep_n_data[:, 12]**2.)
    ep_3_ideal = sqrt(ep_n_data[:, 13]**2. + ep_n_data[:, 14]**2.)
    ep_3_shear = sqrt(ep_n_data[:, 15]**2. + ep_n_data[:, 16]**2.)
    ep_3_full = sqrt(ep_n_data[:, 17]**2. + ep_n_data[:, 18]**2.)

    psi2_org_ideal = compute_psi_n(ep_n_data[:, 1], ep_n_data[:, 2],
                                   ecc_n_data[0, 3], ecc_n_data[0, 4])
    psi2_org_shear = compute_psi_n(ep_n_data[:, 3], ep_n_data[:, 4],
                                   ecc_n_data[0, 3], ecc_n_data[0, 4])
    psi2_org_full = compute_psi_n(ep_n_data[:, 5], ep_n_data[:, 6],
                                  ecc_n_data[0, 3], ecc_n_data[0, 4])
    psi2_ideal = compute_psi_n(ep_n_data[:, 7], ep_n_data[:, 8],
                               ecc_n_data[0, 3], ecc_n_data[0, 4])
    psi2_shear = compute_psi_n(ep_n_data[:, 9], ep_n_data[:, 10],
                               ecc_n_data[0, 3], ecc_n_data[0, 4])
    psi2_full = compute_psi_n(ep_n_data[:, 11], ep_n_data[:, 12],
                              ecc_n_data[0, 3], ecc_n_data[0, 4])
    psi_2_init_list.append([psi2_ideal[0], psi2_shear[0], psi2_full[0]])
    psi_2_final_list.append([psi2_ideal[-1], psi2_shear[-1], psi2_full[-1]])
    psi_2_org_init_list.append([psi2_org_ideal[0], psi2_org_shear[0],
                                psi2_org_full[0]])
    psi_2_org_final_list.append([psi2_org_ideal[-1], psi2_org_shear[-1],
                                 psi2_org_full[-1]])
    psi3_ideal = compute_psi_n(ep_n_data[:, 13], ep_n_data[:, 14],
                               ecc_n_data[0, 5], ecc_n_data[0, 6])
    psi3_shear = compute_psi_n(ep_n_data[:, 15], ep_n_data[:, 16],
                               ecc_n_data[0, 5], ecc_n_data[0, 6])
    psi3_full = compute_psi_n(ep_n_data[:, 17], ep_n_data[:, 18],
                              ecc_n_data[0, 5], ecc_n_data[0, 6])
    psi_3_init_list.append([psi3_ideal[0], psi3_shear[0], psi3_full[0]])
    psi_3_final_list.append([psi3_ideal[-1], psi3_shear[-1], psi3_full[-1]])

    ecc2_dot_ep2_org_ideal = (
        ecc_n_data[:, 3]*ep_n_data[:, 1] + ecc_n_data[:, 4]*ep_n_data[:, 2])
    ecc2_dot_ep2_org_shear = (
        ecc_n_data[:, 3]*ep_n_data[:, 3] + ecc_n_data[:, 4]*ep_n_data[:, 4])
    ecc2_dot_ep2_org_full = (
        ecc_n_data[:, 3]*ep_n_data[:, 5] + ecc_n_data[:, 4]*ep_n_data[:, 6])
    ecc2_dot_ep2_ideal = (
        ecc_n_data[:, 3]*ep_n_data[:, 7] + ecc_n_data[:, 4]*ep_n_data[:, 8])
    ecc2_dot_ep2_shear = (
        ecc_n_data[:, 3]*ep_n_data[:, 9] + ecc_n_data[:, 4]*ep_n_data[:, 10])
    ecc2_dot_ep2_full = (
        ecc_n_data[:, 3]*ep_n_data[:, 11] + ecc_n_data[:, 4]*ep_n_data[:, 12])
    ecc3_dot_ep3_ideal = (
        ecc_n_data[:, 5]*ep_n_data[:, 13] + ecc_n_data[:, 6]*ep_n_data[:, 14])
    ecc3_dot_ep3_shear = (
        ecc_n_data[:, 5]*ep_n_data[:, 15] + ecc_n_data[:, 6]*ep_n_data[:, 16])
    ecc3_dot_ep3_full = (
        ecc_n_data[:, 5]*ep_n_data[:, 17] + ecc_n_data[:, 6]*ep_n_data[:, 18])

    pack_arrays(ep2, ep_2_ideal, ep_2_shear, ep_2_full, ntau, ntau_loc)
    pack_arrays(ep2_org, ep_2_org_ideal, ep_2_org_shear, ep_2_org_full,
                ntau, ntau_loc)
    pack_arrays(ep3, ep_3_ideal, ep_3_shear, ep_3_full, ntau, ntau_loc)
    pack_arrays(psi_2, psi2_ideal, psi2_shear, psi2_full, ntau, ntau_loc)
    pack_arrays(psi_2_org, psi2_org_ideal, psi2_org_shear, psi2_org_full,
                ntau, ntau_loc)
    pack_arrays(psi_3, psi3_ideal, psi3_shear, psi3_full, ntau, ntau_loc)
    pack_arrays(ecc2_dot_ep2, ecc2_dot_ep2_ideal, ecc2_dot_ep2_shear,
                ecc2_dot_ep2_full, ntau, ntau_loc)
    pack_arrays(ecc2_dot_ep2_org, ecc2_dot_ep2_org_ideal,
                ecc2_dot_ep2_org_shear, ecc2_dot_ep2_org_full, ntau, ntau_loc)
    pack_arrays(ecc3_dot_ep3, ecc3_dot_ep3_ideal, ecc3_dot_ep3_shear,
                ecc3_dot_ep3_full, ntau, ntau_loc)

compute_mean_and_stat_error(ep2, nev)
compute_mean_and_stat_error(ep2_org, nev)
compute_mean_and_stat_error(ep3, nev)
compute_mean_and_stat_error(psi_2, nev)
compute_mean_and_stat_error(psi_2_org, nev)
compute_mean_and_stat_error(psi_3, nev)
compute_mean_and_stat_error(ecc2_dot_ep2, nev)
compute_mean_and_stat_error(ecc2_dot_ep2_org, nev)
compute_mean_and_stat_error(ecc3_dot_ep3, nev)

filename = 'momentum_aniso_ep2_evo.dat'
header_text = ("# tau  ep2_ideal  ep2_ideal_err  ep2_shear  ep2_shear_err  "
                + "ep2_full  ep2_full_err")
savetxt(filename, ep2, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'momentum_aniso_ep2_org_evo.dat'
header_text = ("# tau  ep2_ideal  ep2_ideal_err  ep2_shear  ep2_shear_err  "
                + "ep2_full  ep2_full_err")
savetxt(filename, ep2_org, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'momentum_aniso_ep3_evo.dat'
header_text = ("# tau  ep3_ideal  ep3_ideal_err  ep3_shear  ep3_shear_err  "
                + "ep3_full  ep3_full_err")
savetxt(filename, ep3, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'momentum_aniso_psi2_evo.dat'
header_text = ("# tau  psi2_ideal  psi2_ideal_err  "
                + "psi2_shear  psi2_shear_err  psi2_full  psi2_full_err")
savetxt(filename, psi_2, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'momentum_aniso_psi2_org_evo.dat'
header_text = ("# tau  psi2_ideal  psi2_ideal_err  "
                + "psi2_shear  psi2_shear_err  psi2_full  psi2_full_err")
savetxt(filename, psi_2_org, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'momentum_aniso_psi3_evo.dat'
header_text = ("# tau  psi3_ideal  psi3_ideal_err  "
                + "psi3_shear  psi3_shear_err  psi3_full  psi3_full_err")
savetxt(filename, psi_3, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'ecc2_dot_ep2_evo.dat'
header_text = ("# tau  ideal  ideal_err  shear  shear_err  full  full_err")
savetxt(filename, ecc2_dot_ep2, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)

filename = 'ecc2_dot_ep2_org_evo.dat'
header_text = ("# tau  ideal  ideal_err  shear  shear_err  full  full_err")
savetxt(filename, ecc2_dot_ep2_org, fmt="%.8e", delimiter="  ",
        header=header_text)
shutil.move(filename, avg_folder)

filename = 'ecc3_dot_ep3_evo.dat'
header_text = ("# tau  ideal  ideal_err  shear  shear_err  full  full_err")
savetxt(filename, ecc3_dot_ep3, fmt="%.8e", delimiter="  ", header=header_text)
shutil.move(filename, avg_folder)


psi_2_init_list = array(psi_2_init_list)
psi_2_final_list = array(psi_2_final_list)
psi_2_org_init_list = array(psi_2_org_init_list)
psi_2_org_final_list = array(psi_2_org_final_list)
psi_3_init_list = array(psi_3_init_list)
psi_3_final_list = array(psi_3_final_list)


filename_list = ['ideal', 'shear', 'full']
psi_bins = linspace(-pi, pi, 23)
for ifile, filename in enumerate(filename_list):
    hist1, bins1 = histogram(psi_2_init_list[:, ifile], bins=psi_bins)
    hist2, bins2 = histogram(psi_2_final_list[:, ifile], bins=psi_bins)
    width = 0.5*(psi_bins[1] - psi_bins[0])
    center1 = (psi_bins[:-1] + psi_bins[1:])/2 - 0.5*width
    center2 = (psi_bins[:-1] + psi_bins[1:])/2 + 0.5*width
    fig = plt.figure()
    ax = plt.axes([0.12, 0.12, 0.83, 0.83])
    plt.bar(center1, hist1/nev, align='center', width=width, label='initial')
    plt.bar(center2, hist2/nev, align='center', width=width, label='final')
    plt.legend(loc=0, fontsize=18)
    plt.xlabel(r"$2(\Psi_{2p} - \Psi_2)$", fontsize=18)
    plt.ylabel(r"$P(2(\Psi_{2p} - \Psi_2))$", fontsize=18)
    plt.savefig("{0}/{0}_{1}_psi_2_dist.pdf".format(
        results_folder_name, filename), fmt="pdf")
    
    hist1, bins1 = histogram(psi_2_org_init_list[:, ifile], bins=psi_bins)
    hist2, bins2 = histogram(psi_2_org_final_list[:, ifile], bins=psi_bins)
    width = 0.5*(psi_bins[1] - psi_bins[0])
    center1 = (psi_bins[:-1] + psi_bins[1:])/2 - 0.5*width
    center2 = (psi_bins[:-1] + psi_bins[1:])/2 + 0.5*width
    fig = plt.figure()
    ax = plt.axes([0.12, 0.12, 0.83, 0.83])
    plt.bar(center1, hist1/nev, align='center', width=width, label='initial')
    plt.bar(center2, hist2/nev, align='center', width=width, label='final')
    plt.legend(loc=0, fontsize=18)
    plt.xlabel(r"$2(\Psi_{2p} - \Psi_2)$", fontsize=18)
    plt.ylabel(r"$P(2(\Psi_{2p} - \Psi_2))$", fontsize=18)
    plt.savefig("{0}/{0}_{1}_psi_2_org_dist.pdf".format(
        results_folder_name, filename), fmt="pdf")

    hist1, bins1 = histogram(psi_3_init_list[:, ifile], bins=psi_bins)
    hist2, bins2 = histogram(psi_3_final_list[:, ifile], bins=psi_bins)
    width = 0.5*(psi_bins[1] - psi_bins[0])
    center1 = (psi_bins[:-1] + psi_bins[1:])/2 - 0.5*width
    center2 = (psi_bins[:-1] + psi_bins[1:])/2 + 0.5*width
    fig = plt.figure()
    ax = plt.axes([0.14, 0.12, 0.81, 0.83])
    plt.bar(center1, hist1/nev, align='center', width=width, label='initial')
    plt.bar(center2, hist2/nev, align='center', width=width, label='final')
    plt.legend(loc=0, fontsize=18)
    plt.xlabel(r"$3(\Psi_{3p} - \Psi_3)$", fontsize=18)
    plt.ylabel(r"$P(3(\Psi_{3p} - \Psi_3))$", fontsize=18)
    plt.savefig("{0}/{0}_{1}_psi_3_dist.pdf".format(
        results_folder_name, filename), fmt="pdf")
