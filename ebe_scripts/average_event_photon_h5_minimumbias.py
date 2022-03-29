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
centrality_cut_list = Reg_centrality_cut_list + [0., 20., 40., 60., 80.]

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


def check_an_event_is_good(h5_event):
    """This function checks the given event contains all required files"""
    required_files_list = [
        'particle_9999_vndata_eta_-0.5_0.5.dat',
        'particle_9999_vndata_diff_eta_0.1_1.dat',
        'particle_9999_vndata_diff_eta_-1_-0.1.dat',
        'QGP_2to2_total_Spvn_tot_ypTdiff.dat',
        'QGP_AMYcollinear_Spvn_tot_ypTdiff.dat',
        'HG_rho_spectralfun_Spvn_tot_ypTdiff.dat',
        'HG_pipi_bremsstrahlung_Spvn_tot_ypTdiff.dat',
        'HG_omega_Spvn_tot_ypTdiff.dat',
        'HG_2to2_meson_total_Spvn_tot_ypTdiff.dat',
        'photon_total_Spvn.dat',
    ]
    event_file_list = list(h5_event.keys())
    for ifile in required_files_list:
        if ifile not in event_file_list:
            print("event {} is bad, missing {} ...".format(h5_event.name,
                                                           ifile))
            return False
    return True


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
        N_event = data[:, -1]
    else:
        dN_event = data[:, 1]
        N_event = data[:, 1]
    # dN/(pT*dpT*dy)
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    meanpT = sum(dN_interp*pT_inte_array**2.)/sum(dN_interp*pT_inte_array)
    N_interp = exp(interp(pT_inte_array, pT_event, log(N_event+1e-30)))
    N = sum(N_interp)*dpT/0.1
    temp_vn_array = [N, meanpT,]
    for iorder in range(1, n_order):
        if fileType == 0:
            vn_real_event = data[:, 4*iorder]
            vn_imag_event = data[:, 4*iorder+2]
        else:
            vn_real_event = data[:, 2*iorder]
            vn_imag_event = data[:, 2*iorder+1]
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


def calcualte_vn_2_with_gap(vn_data_array, vn_data_array_sub1,
                            vn_data_array_sub2):
    """
        this function computes vn{2} and its stat. err.
        using two subevents with a eta gap
    """
    vn_data_array = array(vn_data_array)
    vn_data_array_sub1 = array(vn_data_array_sub1)
    vn_data_array_sub2 = array(vn_data_array_sub2)
    nev = len(vn_data_array_sub1[:, 0])
    dN = real(vn_data_array[:, 0])
    dN = dN.reshape(len(dN), 1)
    dN1 = real(vn_data_array_sub1[:, 0])
    dN1 = dN1.reshape(len(dN1), 1)
    dN2 = real(vn_data_array_sub2[:, 0])
    dN2 = dN1.reshape(len(dN2), 1)
    Qn_array = dN*vn_data_array[:, 2:]
    Qn_array1 = dN1*vn_data_array_sub1[:, 2:]
    Qn_array2 = dN2*vn_data_array_sub2[:, 2:]
    norder = len(Qn_array[0, :])

    # calcualte observables with Jackknife resampling method
    vnSP_arr = zeros([nev, norder])
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        vnSP_arr[iev, :] = (
            (mean(Qn_array[array_idx, :]*(  conj(Qn_array1[array_idx, :])
                                          + conj(Qn_array2[array_idx, :]))/2.,
                  axis=0))
            /sqrt(mean(dN[array_idx]**2.
                       *real(Qn_array1[array_idx, :]
                             *conj(Qn_array2[array_idx, :])), axis=0)))
    vnSP_mean = real(mean(vnSP_arr, axis=0))
    vnSP_err  = sqrt((nev - 1.)/nev*sum((real(vnSP_arr) - vnSP_mean)**2.,
                                        axis=0))
    return(nan_to_num(vnSP_mean), nan_to_num(vnSP_err))


def calculate_diff_vn_single_event(pT_ref_low, pT_ref_high,
                                   data, data_ref1, data_ref2):
    """
        This function computes pT differential vn{SP} for a single event
        It returns [Qn_pT_arr, Qn_ref_arr]
    """
    npT = 50
    pT_inte_array = linspace(pT_ref_low, pT_ref_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 1]  # photon
    dN_ref1_interp = exp(interp(pT_inte_array, data_ref1[:, 0],
                                log(data_ref1[:, -1] + 1e-30)))
    dN_ref1 = sum(dN_ref1_interp)*dpT/0.1
    dN_ref2_interp = exp(interp(pT_inte_array, data_ref2[:, 0],
                                log(data_ref2[:, -1] + 1e-30)))
    dN_ref2 = sum(dN_ref2_interp)*dpT/0.1
    temp_Qn_pT_array   = [dN_event,]
    temp_Qn_ref1_array = [dN_ref1,]
    temp_Qn_ref2_array = [dN_ref2,]
    for iorder in range(1, n_order):
        vn_real_event = data[:, 2*iorder]      # photon cos
        vn_imag_event = data[:, 2*iorder+1]    # photon sin
        vn_ref_real_interp = interp(pT_inte_array, data_ref1[:, 0],
                                    data_ref1[:, 4*iorder])
        vn_ref_imag_interp = interp(pT_inte_array, data_ref1[:, 0],
                                    data_ref1[:, 4*iorder+2])
        vn_ref_real_inte = (
            sum(vn_ref_real_interp*dN_ref1_interp)/sum(dN_ref1_interp))
        vn_ref_imag_inte = (
            sum(vn_ref_imag_interp*dN_ref1_interp)/sum(dN_ref1_interp))
        Qn_ref1 = dN_ref1*(vn_ref_real_inte + 1j*vn_ref_imag_inte)
        vn_ref_real_interp = interp(pT_inte_array, data_ref2[:, 0],
                                    data_ref2[:, 4*iorder])
        vn_ref_imag_interp = interp(pT_inte_array, data_ref2[:, 0],
                                    data_ref2[:, 4*iorder+2])
        vn_ref_real_inte = (
            sum(vn_ref_real_interp*dN_ref2_interp)/sum(dN_ref2_interp))
        vn_ref_imag_inte = (
            sum(vn_ref_imag_interp*dN_ref2_interp)/sum(dN_ref2_interp))
        Qn_ref2 = dN_ref2*(vn_ref_real_inte + 1j*vn_ref_imag_inte)
        Qn_pt  = dN_event*(vn_real_event + 1j*vn_imag_event)
        temp_Qn_pT_array.append(Qn_pt)
        temp_Qn_ref1_array.append(Qn_ref1)
        temp_Qn_ref2_array.append(Qn_ref2)
    return(temp_Qn_pT_array, temp_Qn_ref1_array, temp_Qn_ref2_array)


def calculate_vn_diff_SP(QnpT_diff, Qnref1, Qnref2,
                         avg_folder, output_filename,
                         pT_array, dN_spectra, dN_spectra_err):
    """
        this funciton calculates the scalar-product vn
        assumption: QnpT is photon, Qnref is charged hadrons
        inputs: QnpT_diff[nev, norder, npT], Qnref[nev, norder]
        output the results to avg_folder/output_filename
        particle spectra result is attached in the output
    """
    QnpT_diff = array(QnpT_diff)
    Qnref1 = array(Qnref1)
    Qnref2 = array(Qnref2)
    nev, norder, npT = QnpT_diff.shape

    vn_diff_SP = []
    Nref1 = real(Qnref1[:, 0])
    Nref2 = real(Qnref2[:, 0])
    N2refPairs = Nref1*Nref2
    NpTPOI = real(QnpT_diff[:, 0, :])
    N2POIPairs = NpTPOI*((Nref1 + Nref2).reshape(nev, 1))
    for iorder in range(1, norder):
        # compute Cn^ref{2}
        QnRef1_tmp = Qnref1[:, iorder]
        QnRef2_tmp = Qnref2[:, iorder]
        n2ref = real(Qnref1[:, iorder]*conj(Qnref2[:, iorder]))

        # compute vn{SP}(pT)
        QnpT_tmp = QnpT_diff[:, iorder, :]
        n2pT = real(  QnpT_tmp*conj(QnRef1_tmp.reshape(nev, 1))
                    + QnpT_tmp*conj(QnRef2_tmp.reshape(nev, 1)))

        # calcualte observables with Jackknife resampling method
        vnSPpT_arr = zeros([nev, npT])
        for iev in range(nev):
            array_idx = [True]*nev
            array_idx[iev] = False
            array_idx = array(array_idx)

            Cn2ref_arr = mean(n2ref[array_idx])/mean(N2refPairs[array_idx])
            vnSPpT_arr[iev, :] = (
                mean(n2pT[array_idx, :], 0)
                /mean(N2POIPairs[array_idx, :], 0)/sqrt(Cn2ref_arr)
            )
        vnSPpT_mean = mean(vnSPpT_arr, 0)
        vnSPpT_err  = sqrt((nev - 1.)/nev
                           *sum((vnSPpT_arr - vnSPpT_mean)**2., 0))
        vn_diff_SP.append(nan_to_num(vnSPpT_mean))
        vn_diff_SP.append(nan_to_num(vnSPpT_err))

    f = open(path.join(avg_folder, output_filename), 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            + "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_array)):
        f.write("%.6e  %.6e  %.6e  "
                % (pT_array[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.6e  %.6e  " % (vn_diff_SP[2*iorder-2][ipT],
                                      vn_diff_SP[2*iorder-1][ipT]))
        f.write("\n")
    f.close()
    return


hf = h5py.File(data_path, "r")
event_list = list(hf.keys())
print("total number of events: {}".format(len(event_list)))

dNdyDict = {}
for ifolder, event_name in enumerate(event_list):
    event_group = hf.get(event_name)
    eventStatus = check_an_event_is_good(event_group)
    file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
    if eventStatus:
        temp_data = event_group.get(file_name)
        temp_data = nan_to_num(temp_data)
        dNdyDict[event_name] = temp_data[0, 1]
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
    print("dNch/deta in [{0:.2f}, {1:.2f}]".format(dN_dy_cut_low,
                                                   dN_dy_cut_high))
    if nev == 0:
        print("Skip ...")
        continue

    print("processing photon ...")

    refFileName0 = 'particle_9999_vndata_diff_eta_-0.5_0.5.dat'
    refFileName1 = 'particle_9999_vndata_diff_eta_0.1_1.dat'
    refFileName2 = 'particle_9999_vndata_diff_eta_-1_-0.1.dat'
    photonFileList = [
        'QGP_2to2_total_Spvn_tot_ypTdiff.dat',
        'QGP_AMYcollinear_Spvn_tot_ypTdiff.dat',
        'HG_rho_spectralfun_Spvn_tot_ypTdiff.dat',
        'HG_pipi_bremsstrahlung_Spvn_tot_ypTdiff.dat',
        'HG_omega_Spvn_tot_ypTdiff.dat',
        'HG_2to2_meson_total_Spvn_tot_ypTdiff.dat',
    ]

    event_group = hf.get(selected_events_list[0])
    temp = nan_to_num(event_group.get("photon_total_Spvn.dat"))
    NPT = temp.shape[0]
    temp = nan_to_num(event_group.get(photonFileList[0]))
    pT_array = temp[:, 1].reshape(-1, NPT)[0, :]
    y_array = temp[:, 0].reshape(-1, NPT)[:, 0]

    for iy, yrap in enumerate(list(y_array)):
        dN_array = []
        vn_phenix_array = []
        vn_phenix_array_ref1 = []; vn_phenix_array_ref2 = []
        vn_alice_array = []
        vn_alice_array_ref1 = []; vn_alice_array_ref2 = []
        QnpT_diff_phenix = []; Qnref1_phenix = []; Qnref2_phenix = []
        QnpT_diff_alice = []; Qnref1_alice = []; Qnref2_alice = []
        for ifolder, event_name in enumerate(selected_events_list):
            event_group    = hf.get(event_name)
            temp_data_ref0 = nan_to_num(event_group.get(refFileName0))
            temp_data_ref1 = nan_to_num(event_group.get(refFileName1))
            temp_data_ref2 = nan_to_num(event_group.get(refFileName2))

            photonRes = []
            for iphoton, photonFilename in enumerate(photonFileList):
                temp_data = nan_to_num(event_group.get(photonFilename))
                nrow, ncol = temp_data.shape

                if iphoton == 0:
                    photonRes = temp_data
                    photonRes[:, 3:] = (temp_data[:, 3:]
                                        *temp_data[:, 2].reshape(nrow, 1))
                else:
                    photonRes[:, 2] += temp_data[:, 2]
                    photonRes[:, 3:] += (temp_data[:, 3:]
                                         *temp_data[:, 2].reshape(nrow, 1))
            photonRes[:, 3:] /= photonRes[:, 2].reshape(nrow, 1)

            photonRes = photonRes[iy*NPT:(iy+1)*NPT, 1:]
            dN_event = photonRes[:, 1]   # dN/(dy pT dpT)

            # record particle spectra
            dN_array.append(dN_event)

            # pT-integrated vn
            # vn with PHENIX pT cut
            temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0, photonRes, 1)
            vn_phenix_array.append(temp_vn_array)
            temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0,
                                                     temp_data_ref1, 0)
            vn_phenix_array_ref1.append(temp_vn_array)
            temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0,
                                                     temp_data_ref2, 0)
            vn_phenix_array_ref2.append(temp_vn_array)

            # vn with ALICE pT cut
            temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0, photonRes, 1)
            vn_alice_array.append(temp_vn_array)
            temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0,
                                                     temp_data_ref1, 0)
            vn_alice_array_ref1.append(temp_vn_array)
            temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0,
                                                     temp_data_ref2, 0)
            vn_alice_array_ref2.append(temp_vn_array)

            # pT-differential vn using scalar-product method
            # vn{SP}(pT) with PHENIX pT cut
            temp_arr = calculate_diff_vn_single_event(0.2, 2.0, photonRes,
                                                      temp_data_ref1,
                                                      temp_data_ref2)
            QnpT_diff_phenix.append(temp_arr[0])
            Qnref1_phenix.append(temp_arr[1])
            Qnref2_phenix.append(temp_arr[2])

            # vn{SP}(pT) with ALICE pT cut
            temp_arr = calculate_diff_vn_single_event(0.20, 3.0, photonRes,
                                                      temp_data_ref1,
                                                      temp_data_ref2)
            QnpT_diff_alice.append(temp_arr[0])
            Qnref1_alice.append(temp_arr[1])
            Qnref2_alice.append(temp_arr[2])

        # now we perform event average
        dN_array = array(dN_array)
        # dN/(2pi dy pT dpT)
        dN_spectra = mean(dN_array, 0)/(2.*pi)
        dN_spectra_err = std(dN_array, 0)/(2.*pi)/sqrt(nev)

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
                vn_phenix_array, vn_phenix_array_ref1, vn_phenix_array_ref2)
        vn_alice_2_gap, vn_alice_2_gap_err = calcualte_vn_2_with_gap(
                vn_alice_array, vn_alice_array_ref1, vn_alice_array_ref2)

        # calcualte vn{SP}(pT)
        output_filename = (
            "photon_differential_observables_PHENIX_y_{}.dat".format(yrap))
        calculate_vn_diff_SP(QnpT_diff_phenix, Qnref1_phenix, Qnref2_phenix,
                             avg_folder, output_filename,
                             pT_array, dN_spectra, dN_spectra_err)
        output_filename = (
            "photon_differential_observables_ALICE_y_{}.dat".format(yrap))
        calculate_vn_diff_SP(QnpT_diff_alice, Qnref1_alice, Qnref2_alice,
                             avg_folder, output_filename,
                             pT_array, dN_spectra, dN_spectra_err)

        ######################################################################
        # finally, output all the results
        ######################################################################

        output_filename = "photon_integrated_observables_y_{}.dat".format(yrap)
        f = open(path.join(avg_folder, output_filename), 'w')
        f.write("dN/dy(phenix)= %.5e +/- %.5e\n" % (dNdy_avg_phenix,
                                                    dNdy_avg_phenix_err))
        f.write("dN/dy(ALICE)= %.5e +/- %.5e\n" % (dNdy_avg_alice,
                                                   dNdy_avg_alice_err))
        f.write("<pT>(phenix)= %.5e +/- %.5e\n" % (meanpT_phenix,
                                                   meanpT_phenix_err))
        f.write("<pT>(ALICE)= %.5e +/- %.5e\n" % (meanpT_alice,
                                                  meanpT_alice_err))
        for iorder in range(1, n_order):
            f.write("v_%d{2}(phenix)= %.5e +/- %.5e\n"
                    % (iorder, vn_phenix_2_gap[iorder-1],
                       vn_phenix_2_gap_err[iorder-1]))
            f.write("v_%d{2}(ALICE)= %.5e +/- %.5e\n"
                    % (iorder, vn_alice_2_gap[iorder-1],
                       vn_alice_2_gap_err[iorder-1]))
        f.close()
print("Analysis is done.")
