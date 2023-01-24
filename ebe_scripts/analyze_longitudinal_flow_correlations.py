#! /usr/bin/env python3
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
from scipy.optimize import curve_fit

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
centralityCutList = Reg_centrality_cut_list

CentralityFlag = 1   # 1: sort dNch/deta and cut centrality

RapidityTrigger = 0  # 0: mid-rapidity [-0.5, 0.5]
                     # 1: PHENIX BBC trigger [-3.9, -3.1]
                     # 2: ALICE V0A trigger [-5.1, -2.8]
                     # 3: ATLAS forward trigger [-4.9, -3.1]

RapTrigLabel = "CL1"
if RapidityTrigger == 1:
    RapTrigLabel = "BBC"
elif RapidityTrigger == 2:
    RapTrigLabel = "V0A"
elif RapidityTrigger == 3:
    RapTrigLabel = "ATLASForward"

try:
    data_path = path.abspath(argv[1])
    data_name = data_path.split("/")[-1]
    resultsFolderName = data_name.split(".h5")[0]
    avg_folder_header = path.join(
        path.abspath(argv[2]), "{}_{}".format(resultsFolderName, RapTrigLabel))
    print("output folder: %s" % avg_folder_header)
    if path.isdir(avg_folder_header):
        print("folder %s already exists!" % avg_folder_header)
        var = input("do you want to delete it? [y/N]")
        if 'y' in var.lower():
            shutil.rmtree(avg_folder_header)
            mkdir(avg_folder_header)
        else:
            print("Continue analysis in {} ...".format(avg_folder_header))
    else:
        mkdir(avg_folder_header)
except IndexError:
    print("Usage: {} <path_to_hdf5_dataFile> <results_folder>".format(argv[0]))
    exit(1)

particle_list = ['9999']
particle_name_list = ['charged_hadron']

n_order = 7


def LinearFunc(x, a, b):
    return a * x + b


def calculate_vn_eta(eta_array, dN_array, vn_array, eta_min, eta_max):
    """
        This function computes vn(eta).
        eta_min and eta_max specify the rapidity range of reference flow vector
    """
    nev, neta   = dN_array.shape
    dN_array    = dN_array.reshape((nev, 1, neta))
    idx         = (eta_array > eta_min) & (eta_array < eta_max)
    vn_ref      = (sum(dN_array[:, :, idx]*vn_array[:, :, idx], axis=2)
                   /(sum(dN_array[:, :, idx], axis=2) + 1e-15))
    vnshape     = vn_ref.shape
    nvn         = vnshape[1]
    vn_ref      = vn_ref.reshape((vnshape[0], vnshape[1], 1))
    vn_SP_ev    = real(vn_array*conj(vn_ref))
    vn_SP_array = zeros([nev, nvn, neta])
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        vn_den         = mean((absolute(vn_ref[array_idx, :, :]))**2., axis=0)
        vn_SP          = (mean(vn_SP_ev[array_idx, :, :], axis=0)
                          /(sqrt(vn_den) + 1e-16))
        vn_SP_array[iev, :, :] = vn_SP
    vn_SP_mean = mean(vn_SP_array, axis=0)
    vn_SP_err  = sqrt((nev - 1.)/nev*sum((vn_SP_array - vn_SP_mean)**2., axis=0))
    return([vn_SP_mean, vn_SP_err])


def calculate_rn_eta(eta_array, eta_min, eta_max, dN_array, vn_array,
                     outputFileName):
    """
        This function computes the longitudinal factorization breaking ratios
        for all n passed from vn_array
            eta, rn(eta), r_nn(eta)
        the reference flow angles are defined in [eta_min, eta_max]
        and [-eta_max, -eta_min]
    """
    nev, neta = dN_array.shape
    dN_array = dN_array.reshape((nev, 1, neta))
    Qn_array = vn_array
    nQn = Qn_array.shape[1]

    # calculate the reference flow vector for every event
    eta_b_min    = abs(eta_min)
    eta_b_max    = abs(eta_max)
    eta_ref1_tmp = linspace(eta_b_min, eta_b_max, 16)
    eta_ref2_tmp = linspace(-eta_b_max, -eta_b_min, 16)
    Qn_ref1      = []
    Qn_ref2      = []
    for iev in range(nev):
        dN1_interp = interp(eta_ref1_tmp, eta_array, dN_array[iev, 0, :])
        dN2_interp = interp(eta_ref2_tmp, eta_array, dN_array[iev, 0, :])
        Qn_ref1_vec = []
        Qn_ref2_vec = []
        for iorder in range(nQn):
            Qn1_interp = interp(eta_ref1_tmp, eta_array,
                                Qn_array[iev, iorder, :])
            Qn2_interp = interp(eta_ref2_tmp, eta_array,
                                Qn_array[iev, iorder, :])
            Qn_ref1_vec.append(sum(dN1_interp*Qn1_interp)
                               /(sum(dN1_interp) + 1e-15))
            Qn_ref2_vec.append(sum(dN2_interp*Qn2_interp)
                               /(sum(dN2_interp) + 1e-15))
        Qn_ref1.append(Qn_ref1_vec)
        Qn_ref2.append(Qn_ref2_vec)
    Qn_ref1 = array(Qn_ref1).reshape((nev, nQn, 1))
    Qn_ref2 = array(Qn_ref2).reshape((nev, nQn, 1))

    rn_num    = real(Qn_array[:, :, ::-1]*conj(Qn_ref1))
    rn_den    = real(Qn_array*conj(Qn_ref1))
    rnn_num   = real((Qn_ref2*conj(Qn_array))
                     *(Qn_array[:, :, ::-1]*conj(Qn_ref1)))
    rnn_den   = real((Qn_ref2*conj(Qn_array[:, :, ::-1]))
                     *(Qn_array*conj(Qn_ref1)))

    # compute the error using jack-knife
    rn_array  = zeros([nev, nQn, neta])
    rnn_array = zeros([nev, nQn, neta])
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        rn_ev          = (mean(rn_num[array_idx], axis=0)
                          /(mean(rn_den[array_idx], axis=0) + 1e-15))
        rnn_ev         = (mean(rnn_num[array_idx], axis=0)
                          /(mean(rnn_den[array_idx], axis=0) + 1e-15))
        rn_array[iev, :, :]  = rn_ev
        rnn_array[iev, :, :] = rnn_ev
    rn_mean  = mean(rn_array, axis=0)
    rn_err   = sqrt((nev - 1.)/nev*sum((rn_array - rn_mean)**2., axis=0))
    rnn_mean = mean(rnn_array, axis=0)
    rnn_err  = sqrt((nev - 1.)/nev*sum((rnn_array - rnn_mean)**2., axis=0))

    # perform a linear fit for r_2, r_3, and r_4
    rn_slope = []
    rnn_slope = []
    idx = abs(eta_array) < 1.
    for iorder in range(1, 4):
        popt, pcov = curve_fit(LinearFunc, eta_array[idx],
                               rn_mean[iorder, idx],
                               sigma = rn_err[iorder, idx],
                               method = "dogbox")
        rn_slope.append([popt[0], sqrt(pcov[0, 0])])
        popt, pcov = curve_fit(LinearFunc, eta_array[idx],
                               rnn_mean[iorder, idx],
                               sigma = rnn_err[iorder, idx],
                               method = "dogbox")
        rnn_slope.append([popt[0], sqrt(pcov[0, 0])])

    f = open(outputFileName, 'w')
    f.write("#eta  rn(eta)  rn_err(eta)  rnn(eta)  rnn_err(eta)\n")
    for ieta in range(len(eta_array)-1):
        f.write("%.10e  " % eta_array[ieta])
        for iorder in range(nQn):
            f.write("%.10e  %.10e  %.10e  %.10e  "
                    % (rn_mean[iorder, ieta], rn_err[iorder, ieta],
                       rnn_mean[iorder, ieta], rnn_err[iorder, ieta]))
        f.write("\n")
    f.close()
    return rn_slope, rnn_slope



###############################################################################
### analysis starts from here ...
###############################################################################

hf = h5py.File(data_path, "r")
event_list = list(hf.keys())
print("total number of events: {}".format(len(event_list)))

dNdyDict = {}
dNdyList = []
for ifolder, event_name in enumerate(event_list):
    file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
    if RapidityTrigger == 1:      # PHENIX BBC Trigger
        file_name = "particle_9999_vndata_eta_-3.9_-3.1.dat"
    elif RapidityTrigger == 2:    # ALICE V0A Trigger
        file_name = "particle_9999_vndata_eta_-5.1_-2.8.dat"
    elif RapidityTrigger == 3:    # ATLAS forward Trigger
        file_name = "particle_9999_vndata_eta_-4.9_-3.1.dat"
    try:
        event_group = hf.get(event_name)
        temp_data   = event_group.get(file_name)
        temp_data   = nan_to_num(temp_data)
        dNdyDict[event_name] = temp_data[0, 1]
    except:
        continue
dNdyList = -sort(-array(list(dNdyDict.values())))
print("Number of good events: {}".format(len(dNdyList)))


rnSlopeFile = open(path.join(avg_folder_header, "rnSlope.dat"), "w")
rnSlopeFile.write("#cen  rn_slope  rn_slope_err (n = 2 - 4)\n")
rnnSlopeFile = open(path.join(avg_folder_header, "rnnSlope.dat"), "w")
rnnSlopeFile.write("#cen  rnn_slope  rnn_slope_err (n = 2 - 4)\n")
ecc_rnSlopeFile = open(path.join(avg_folder_header, "ecc_rnSlope.dat"), "w")
ecc_rnSlopeFile.write("#cen  rn_slope  rn_slope_err (n = 2 - 4)\n")
ecc_rnnSlopeFile = open(path.join(avg_folder_header, "ecc_rnnSlope.dat"), "w")
ecc_rnnSlopeFile.write("#cen  rnn_slope  rnn_slope_err (n = 2 - 4)\n")
for icen in range(len(centralityCutList) - 1):
    if centralityCutList[icen+1] < centralityCutList[icen]: continue
    avg_folder = path.join(
        avg_folder_header, "{0:02.0f}-{1:02.0f}".format(
            centralityCutList[icen], centralityCutList[icen+1])
    )

    if path.isdir(avg_folder):
        print("{} already exists, skipped ...".format(avg_folder))
        continue
    else:
        mkdir(avg_folder)

    selected_events_list = []
    dN_dy_cut_high = dNdyList[
        int(len(dNdyList)*centralityCutList[icen]/100.)
    ]
    dN_dy_cut_low  = dNdyList[
        min(len(dNdyList)-1,
            int(len(dNdyList)*centralityCutList[icen+1]/100.))
    ]
    for event_name in dNdyDict.keys():
        if (dNdyDict[event_name] > dN_dy_cut_low
            and dNdyDict[event_name] <= dN_dy_cut_high):
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    print("analysis {}%-{}% nev = {}...".format(
            centralityCutList[icen], centralityCutList[icen+1], nev))
    if nev == 0:
        print("Skip ...")
        continue

    for ipart, particle_id in enumerate(particle_list):
        print("processing %s ..." % particle_name_list[ipart])

        # particle rapidity distribution
        if particle_id == '9999':
            file_name = 'particle_%s_dNdeta_pT_0.2_3.dat' % particle_id
        else:
            file_name = 'particle_%s_dNdy_pT_0.2_3.dat' % particle_id

        eta_array = []
        dN_array = []
        totalN_array = []
        vn_array = []
        ecc_eta_arr = []
        dEdetas_array = []
        eccn_array = []
        for ifolder, event_name in enumerate(selected_events_list):
            event_group = hf.get(event_name)
            eccFileList = []
            for ifile in event_group.keys():
                if "eccentricities_evo_ed" in ifile:
                    eccFileList.append(ifile)
            ecc_data = nan_to_num(event_group.get(eccFileList[0]))
            temp_data = nan_to_num(event_group.get(file_name))

            if ifolder == 0:
                ecc_eta_arr = ecc_data[:, 0]
            dEdetas_array.append(ecc_data[:, 1])
            temp_eccn_array = []
            for iorder in range(1, 6):
                eccn = ecc_data[:, 2*iorder] + 1j*ecc_data[:, 2*iorder+1]
                temp_eccn_array.append(eccn)
            eccn_array.append(temp_eccn_array)

            eta_array.append(temp_data[:, 0])
            dN_array.append(temp_data[:, 1])
            totalN_array.append(temp_data[:, -1])
            temp_vn_array = []
            for iorder in range(1, n_order):
                vn_real = temp_data[:, 2*iorder+1]
                vn_imag = temp_data[:, 2*iorder+2]
                vn = vn_real + 1j*vn_imag
                temp_vn_array.append(vn)
            vn_array.append(temp_vn_array)

        eta_array = array(eta_array)
        dN_array = array(dN_array)
        totalN_array = array(totalN_array)
        vn_array = array(vn_array)
        dEdetas_array = array(dEdetas_array)
        eccn_array = array(eccn_array)

        eta_point = mean(eta_array, 0)
        dNdeta = mean(dN_array, 0)
        dNdeta_err = std(dN_array, 0)/sqrt(nev)
        if RapidityTrigger == 0:
            vn_SP_eta, vn_SP_eta_err = calculate_vn_eta(
                                eta_point, dN_array, vn_array, -3.0, 3.0)
        elif RapidityTrigger == 1:
            vn_SP_eta, vn_SP_eta_err = calculate_vn_eta(
                                eta_point, dN_array, vn_array, -3.9, -3.1)
        elif RapidityTrigger == 2:
            vn_SP_eta, vn_SP_eta_err = calculate_vn_eta(
                                eta_point, dN_array, vn_array, -5.1, -2.8)
        elif RapidityTrigger == 3:
            vn_SP_eta, vn_SP_eta_err = calculate_vn_eta(
                                eta_point, dN_array, vn_array, -4.9, -3.1)
        vn_SP_eta_mid, vn_SP_eta_mid_err = calculate_vn_eta(
                                eta_point, dN_array, vn_array, -0.8, 0.8)

        output_filename = "{}_rapidity_distribution.dat".format(
                                                    particle_name_list[ipart])
        f = open(path.join(avg_folder, output_filename), 'w')
        if(particle_id == '9999'):
            f.write("#eta  dN/deta  dN/deta_err  vn{2}(eta)  vn{2}(eta)_err"
                    + "  Re{vn}(eta) Re{vn}(eta)_err\n")
        else:
            f.write("#y  dN/dy  dN/dy_err  vn{2}(y)  vn{2}(y)_err  "
                    + "Re{vn}(y)  Re{vn}(y)_err\n")
        for ieta in range(len(eta_point)):
            f.write("%.10e  %.10e  %.10e  "
                    % (eta_point[ieta], dNdeta[ieta], dNdeta_err[ieta]))
            for iorder in range(1, n_order):
                f.write("%.10e  %.10e  %.10e  %.10e  "
                        % (vn_SP_eta[iorder-1, ieta],
                           vn_SP_eta_err[iorder-1, ieta],
                           vn_SP_eta_mid[iorder-1, ieta],
                           vn_SP_eta_mid_err[iorder-1, ieta]))
            f.write("\n")
        f.close()

        # calculate the longitudinal flow decorrelation with ATLAS cut
        etaMin = 4.0; etaMax = 4.9
        output_filename = path.join(avg_folder, "ecc_rn_eta.dat")
        ecc_rnSlope, ecc_rnnSlope = calculate_rn_eta(
                ecc_eta_arr, etaMin, etaMax, dEdetas_array, eccn_array,
                output_filename)
        output_filename = path.join(avg_folder, "{}_rn_eta.dat".format(
                                                    particle_name_list[ipart]))
        rnSlope, rnnSlope = calculate_rn_eta(
                eta_point, etaMin, etaMax, dN_array, vn_array, output_filename)
        centralityCenter = (centralityCutList[icen]
                            + centralityCutList[icen+1])/2.
        ecc_rnSlopeFile.write("%.2f  " % centralityCenter)
        for ir, ir_err in ecc_rnSlope:
            ecc_rnSlopeFile.write("%.6e  %.6e  " % (ir, ir_err))
        ecc_rnSlopeFile.write("\n")
        ecc_rnnSlopeFile.write("%.2f  " % centralityCenter)
        for ir, ir_err in ecc_rnnSlope:
            ecc_rnnSlopeFile.write("%.6e  %.6e  " % (ir, ir_err))
        ecc_rnnSlopeFile.write("\n")
        rnSlopeFile.write("%.2f  " % centralityCenter)
        for ir, ir_err in rnSlope:
            rnSlopeFile.write("%.6e  %.6e  " % (ir, ir_err))
        rnSlopeFile.write("\n")
        rnnSlopeFile.write("%.2f  " % centralityCenter)
        for ir, ir_err in rnnSlope:
            rnnSlopeFile.write("%.6e  %.6e  " % (ir, ir_err))
        rnnSlopeFile.write("\n")

rnSlopeFile.close()
rnnSlopeFile.close()
ecc_rnSlopeFile.close()
ecc_rnnSlopeFile.close()
print("Analysis is done.")

