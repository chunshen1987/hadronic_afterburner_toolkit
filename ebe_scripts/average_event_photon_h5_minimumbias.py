#!/usr/bin/env python3
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

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
PHOBOS_cen_list = [0., 6., 15., 25., 35., 45., 55.]  # PHOBOS AuAu 200
SPS_cen_list    = [5., 12.5, 23.5, 33.5, 43.5]       # SPS PbPb
PHENIX_cen_list = [20., 40., 60., 88.]               # PHENIX dAu
STAR_cen_list   = [0., 10., 40., 80]                 # STAR v1
centrality_cut_list = Reg_centrality_cut_list

n_order = 7

try:
    data_path = path.abspath(argv[1])
    data_name = data_path.split("/")[-1]
    results_folder_name = data_name.split(".h5")[0]
    avg_folder_header = path.join(path.abspath(argv[2]),
                                  results_folder_name)
    print("output folder: %s" % avg_folder_header)
    if(path.isdir(avg_folder_header)):
        print("folder %s already exists!" % avg_folder_header)
        var = input("do you want to delete it? [y/N]")
        if 'y' in var:
            shutil.rmtree(avg_folder_header)
        else:
            print("please choose another folder path~")
            exit(0)
    mkdir(avg_folder_header)
except IndexError:
    print("Usage: {} database.h5 results_folder".format(argv[0]))
    exit(1)


def calculate_meanpT_inte_vn(pT_low, pT_high, data, fileType):
    """
        this function calculates the dN/dy, <pT>, pT-integrated vn in a
        given pT range (pT_low, pT_high) for every event in the data
        fileType == 0: hadrons
        fileType == 1: photons
    """
    npT = 50
    pT_inte_array = linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    pT_event = data[:, 0]
    if fileType == 0:
        dN_event = data[:, 2]
    else:
        dN_event = data[:, 1]
    # dN/(2pi*pT*dpT*dy)
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    N = sum(dN_interp*pT_inte_array)*dpT*2.*pi
    meanpT = sum(dN_interp*pT_inte_array**2.)/sum(dN_interp*pT_inte_array)
    temp_vn_array = [N, meanpT,]
    for iorder in range(1, n_order):
        if fileType == 0:
            vn_real_event = data[:, 4*iorder]
            vn_imag_event = data[:, 4*iorder+2]
        else:
            vn_real_event = data[:, 3*iorder-1]
            vn_imag_event = data[:, 3*iorder]
        vn_real_interp = interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = interp(pT_inte_array, pT_event, vn_imag_event)
        vn_real_inte = (
            sum(vn_real_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_imag_inte = (
            sum(vn_imag_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_inte = vn_real_inte + 1j*vn_imag_inte
        temp_vn_array.append(vn_inte)
    return(temp_vn_array)


def calcualte_vn_2_with_gap(vn_data_array_sub1, vn_data_array_sub2):
    """
        this function computes vn{2} and its stat. err.
        using two subevents with a eta gap
    """
    vn_data_array_sub1 = array(vn_data_array_sub1)
    vn_data_array_sub2 = array(vn_data_array_sub2)
    nev = len(vn_data_array_sub1[:, 0])
    dN1 = real(vn_data_array_sub1[:, 0])
    dN1 = dN1.reshape(len(dN1), 1)
    dN2 = real(vn_data_array_sub2[:, 0])
    dN2 = dN1.reshape(len(dN2), 1)
    Qn_array1 = dN1*vn_data_array_sub1[:, 2:]
    Qn_array2 = dN2*vn_data_array_sub2[:, 2:]

    num = sqrt(mean(real(Qn_array1*conj(Qn_array2)), axis=0))
    num_err = std(real(Qn_array1*conj(Qn_array2)), axis=0)/sqrt(nev)/(2.*num)
    denorm = sqrt(mean(dN1*dN2))
    denorm_err = std(dN1*dN2)/sqrt(nev)/(2.*denorm)
    vn_2 = num/denorm
    vn_2_err = sqrt((num_err/denorm)**2. + (num*denorm_err/denorm**2.)**2.)
    return(nan_to_num(vn_2), nan_to_num(vn_2_err))


def calculate_diff_vn_single_event(pT_ref_low, pT_ref_high, data, data_ref):
    """
        This function computes pT differential vn{4} for a single event
        It returns [Qn_pT_arr, Qn_ref_arr]
    """
    npT = 50
    pT_inte_array = linspace(pT_ref_low, pT_ref_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 1]  # photon
    dN_ref_event = data_ref[:, 2]
    pT_ref_event = data_ref[:, 0]
    dN_ref_interp = exp(interp(pT_inte_array, pT_ref_event,
                               log(dN_ref_event + 1e-30)))
    dN_ref = sum(dN_ref_interp*pT_inte_array)*dpT*2.*pi
    temp_Qn_pT_array = [dN_event,]
    temp_Qn_ref_array = [dN_ref,]
    for iorder in range(1, n_order):
        vn_real_event = data[:, 3*iorder-1]  # photon
        vn_imag_event = data[:, 3*iorder]    # photon
        vn_ref_real_event = data_ref[:, 4*iorder]
        vn_ref_imag_event = data_ref[:, 4*iorder+2]
        vn_ref_real_interp = interp(pT_inte_array, pT_ref_event,
                                    vn_ref_real_event)
        vn_ref_imag_interp = interp(pT_inte_array, pT_ref_event,
                                    vn_ref_imag_event)
        vn_ref_real_inte = (
            sum(vn_ref_real_interp*dN_ref_interp)/sum(dN_ref_interp))
        vn_ref_imag_inte = (
            sum(vn_ref_imag_interp*dN_ref_interp)/sum(dN_ref_interp))
        Qn_ref = dN_ref*(vn_ref_real_inte + 1j*vn_ref_imag_inte)
        Qn_pt  = dN_event*(vn_real_event + 1j*vn_imag_event)
        temp_Qn_pT_array.append(Qn_pt)
        temp_Qn_ref_array.append(Qn_ref)
    return(temp_Qn_pT_array, temp_Qn_ref_array)


def calculate_vn_diff_SP(QnpT_diff, Qnref):
    """
        this funciton calculates the scalar-product vn
        assumption: QnpT is photon, Qnref is charged hadrons
        inputs: QnpT_diff[nev, norder, npT], Qnref[nev, norder]
        return: [vn{SP}(pT), vn{SP}(pT)_err]
    """
    QnpT_diff = array(QnpT_diff)
    Qnref = array(Qnref)
    nev, norder, npT = QnpT_diff.shape

    vn_diff_SP = []
    Nref = real(Qnref[:, 0])
    N2refPairs = Nref*(Nref - 1.)
    for iorder in range(1, norder):
        # compute Cn^ref{2}
        QnRef_tmp = Qnref[:, iorder]
        n2ref = abs(QnRef_tmp)**2. - Nref

        # calcualte observables with Jackknife resampling method
        Cn2ref_arr = zeros(nev)
        for iev in range(nev):
            array_idx = [True]*nev
            array_idx[iev] = False
            array_idx = array(array_idx)

            Cn2ref_arr[iev] = (
                    mean(n2ref[array_idx])/mean(N2refPairs[array_idx]))
        Cn2ref_mean = mean(Cn2ref_arr)
        Cn2ref_err  = sqrt((nev - 1.)/nev*sum((Cn2ref_arr - Cn2ref_mean)**2.))

        # compute vn{SP}(pT)
        NpTPOI = real(QnpT_diff[:, 0, :])
        QnpT_tmp = QnpT_diff[:, iorder, :]
        Nref = Nref.reshape(nev, 1)
        QnRef_tmp = QnRef_tmp.reshape(nev, 1)
        N2POIPairs = NpTPOI*Nref
        n2pT = real(QnpT_tmp*conj(QnRef_tmp))

        # calcualte observables with Jackknife resampling method
        vnSPpT_arr = zeros([nev, npT])
        for iev in range(nev):
            array_idx = [True]*nev
            array_idx[iev] = False
            array_idx = array(array_idx)

            vnSPpT_arr[iev, :] = (mean(n2pT[array_idx], 0)
                    /mean(N2POIPairs[array_idx], 0)/sqrt(Cn2ref_arr[iev]))
        vnSPpT_mean = mean(vnSPpT_arr, 0)
        vnSPpT_err  = sqrt((nev - 1.)/nev
                           *sum((vnSPpT_arr - vnSPpT_mean)**2., 0))
        vn_diff_SP.append(vnSPpT_mean)
        vn_diff_SP.append(vnSPpT_err)
    return vn_diff_SP


hf = h5py.File(data_path, "r")
event_list = list(hf.keys())
print("total number of events: {}".format(len(event_list)))

dNdyDict = {}
for ifolder, event_name in enumerate(event_list):
    file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
    try:
        event_group = hf.get(event_name)
        temp_data   = event_group.get(file_name)
        temp_data   = nan_to_num(temp_data)
        dNdyDict[event_name] = temp_data[0, 1]
    except:
        continue
dNdyList = -sort(-array(list(dNdyDict.values())))
print("Number of good events: {}".format(len(dNdyList)))


for icen in range(len(centrality_cut_list) - 1):
    if centrality_cut_list[icen+1] < centrality_cut_list[icen]: continue
    avg_folder = path.join(
        avg_folder_header, "{0:02.0f}-{1:02.0f}".format(
            centrality_cut_list[icen], centrality_cut_list[icen+1])
    )
    mkdir(avg_folder)

    selected_events_list = []
    dN_dy_cut_high = dNdyList[
        int(len(dNdyList)*centrality_cut_list[icen]/100.)
    ]
    dN_dy_cut_low  = dNdyList[
        min(len(dNdyList) - 1,
            int(len(dNdyList)*centrality_cut_list[icen+1]/100.))
    ]
    for event_name in dNdyDict.keys():
        if (dNdyDict[event_name] > dN_dy_cut_low
            and dNdyDict[event_name] <= dN_dy_cut_high):
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    print("analysis {}%-{}% nev = {}...".format(
            centrality_cut_list[icen], centrality_cut_list[icen+1], nev))
    if nev == 0:
        print("Skip ...")
        continue

    print("processing photon ...")

    refFileName = 'particle_9999_vndata_diff_eta_0.5_2.dat'
    photonFilename = 'photon_total_Spvn.dat'

    pT_array = []
    dN_array = []
    vn_phenix_array = []
    vn_phenix_array_ref = []
    vn_alice_array = []
    vn_alice_array_ref = []
    QnpT_diff_phenix = []; Qnref_phenix = []
    QnpT_diff_alice = []; Qnref_alice = []
    for ifolder, event_name in enumerate(selected_events_list):
        event_group   = hf.get(event_name)
        temp_data     = nan_to_num(event_group.get(photonFilename))
        temp_data_ref = nan_to_num(event_group.get(refFileName))

        dN_event = temp_data[:, 1]  # dN/(2pi dy pT dpT)
        pT_event = temp_data[:, 0]

        # record particle spectra
        pT_array.append(pT_event)
        dN_array.append(dN_event)

        # pT-integrated vn
        # vn with PHENIX pT cut
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0, temp_data, 1)
        vn_phenix_array.append(temp_vn_array)
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0, temp_data_ref, 0)
        vn_phenix_array_ref.append(temp_vn_array)

        # vn with ALICE pT cut
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0, temp_data, 1)
        vn_alice_array.append(temp_vn_array)
        temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0, temp_data_ref, 0)
        vn_alice_array_ref.append(temp_vn_array)

        # pT-differential vn using scalar-product method
        # vn{SP}(pT) with PHENIX pT cut
        temp_arr = calculate_diff_vn_single_event(0.2, 2.0, temp_data,
                                                  temp_data_ref)
        QnpT_diff_phenix.append(temp_arr[0])
        Qnref_phenix.append(temp_arr[1])

        # vn{SP}(pT) with ALICE pT cut
        temp_arr = calculate_diff_vn_single_event(0.20, 3.0, temp_data,
                                                  temp_data_ref)
        QnpT_diff_alice.append(temp_arr[0])
        Qnref_alice.append(temp_arr[1])

    # now we perform event average
    dN_array = array(dN_array)
    pT_array = array(pT_array)

    n_pT = len(pT_array[0, :])
    pT_spectra = zeros([n_pT])
    for ipT in range(len(pT_array[0, :])):
        dN_temp = sum(dN_array[:, ipT]*pT_array[:, ipT])
        if dN_temp > 0:
            pT_spectra[ipT] = (
                    sum(pT_array[:, ipT]**2.*dN_array[:, ipT])/dN_temp)
        else:
            pT_spectra[ipT] = mean(pT_array[:, ipT])
    # dN/(2pi dy pT dpT)
    dN_spectra = mean(pT_array*dN_array, 0)/pT_spectra
    dN_spectra_err = std(pT_array*dN_array, 0)/pT_spectra/sqrt(nev)

    # calcualte dN/dy and <pT>
    vn_phenix_array = array(vn_phenix_array)
    vn_alice_array = array(vn_alice_array)
    dNdy_avg_phenix     = real(mean(vn_phenix_array[:, 0]))
    dNdy_avg_phenix_err = real(std(vn_phenix_array[:, 0]))/sqrt(nev)
    meanpT_phenix       = real(mean(vn_phenix_array[:, 1]))
    meanpT_phenix_err   = real(std(vn_phenix_array[:, 1]))/sqrt(nev)
    dNdy_avg_alice     = real(mean(vn_alice_array[:, 0]))
    dNdy_avg_alice_err = real(std(vn_alice_array[:, 0]))/sqrt(nev)
    meanpT_alice       = real(mean(vn_alice_array[:, 1]))
    meanpT_alice_err   = real(std(vn_alice_array[:, 1]))/sqrt(nev)

    # calcualte vn{2}
    vn_phenix_2_gap, vn_phenix_2_gap_err = calcualte_vn_2_with_gap(
                                vn_phenix_array, vn_phenix_array_ref)
    vn_alice_2_gap, vn_alice_2_gap_err = calcualte_vn_2_with_gap(
                                vn_alice_array, vn_alice_array_ref)

    # calcualte vn{SP}(pT)
    vn_diff_SP_phenix = calculate_vn_diff_SP(QnpT_diff_phenix,
                                             Qnref_phenix)
    vn_diff_SP_alice = calculate_vn_diff_SP(QnpT_diff_alice, Qnref_alice)

    ######################################################################
    # finally, output all the results
    ######################################################################

    output_filename = "photon_integrated_observables.dat"
    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("dN/dy(phenix)= %.5e +/- %.5e\n" % (dNdy_avg_phenix,
                                                dNdy_avg_phenix_err))
    f.write("dN/dy(ALICE)= %.5e +/- %.5e\n" % (dNdy_avg_alice,
                                               dNdy_avg_alice_err))
    f.write("<pT>(phenix)= %.5e +/- %.5e\n" % (meanpT_phenix,
                                               meanpT_phenix_err))
    f.write("<pT>(ALICE)= %.5e +/- %.5e\n" % (meanpT_alice, meanpT_alice_err))
    for iorder in range(1, n_order):
        f.write("v_%d{2}(phenix)= %.5e +/- %.5e\n"
                % (iorder, vn_phenix_2_gap[iorder-1],
                   vn_phenix_2_gap_err[iorder-1]))
        f.write("v_%d{2}(ALICE)= %.5e +/- %.5e\n"
                % (iorder, vn_alice_2_gap[iorder-1],
                   vn_alice_2_gap_err[iorder-1]))
    f.close()

    output_filename = "photon_differential_observables_PHENIX.dat"
    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            + "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_phenix[2*iorder-2][ipT],
                                        vn_diff_SP_phenix[2*iorder-1][ipT]))
        f.write("\n")
    f.close()

    output_filename = "photon_differential_observables_ALICE.dat"
    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            + "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_alice[2*iorder-2][ipT],
                                        vn_diff_SP_alice[2*iorder-1][ipT]))
        f.write("\n")
    f.close()

print("Analysis is done.")
