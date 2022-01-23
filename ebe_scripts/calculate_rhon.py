#! /usr/bin/env python3
"""
     This script get event-by-event vn and initial ecc_n
     v_n is analyzed up to n = 6
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
from numpy import *
import h5py
import shutil

Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
PHOBOS_cen_list = [0., 6., 15., 25., 35., 45., 55.]  # PHOBOS AuAu 200
SPS_cen_list    = [5., 12.5, 23.5, 33.5, 43.5]       # SPS PbPb
PHENIX_cen_list = [20., 40., 60., 88.]               # PHENIX dAu
STAR_cen_list   = [0., 10., 40., 80]                 # STAR v1
#centrality_cut_list = (Reg_centrality_cut_list + PHOBOS_cen_list
#                       + SPS_cen_list + PHENIX_cen_list + STAR_cen_list)
centrality_cut_list = Reg_centrality_cut_list
#centrality_cut_list = linspace(0., 100., 21)

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
    print("Usage: {} working_folder results_folder".format(argv[0]))
    exit(1)

particle_list = ['9999', '211', '321', '2212', '-211', '-321', '-2212',
                 '3122', '-3122', '3312', '-3312', '3334', '-3334',
                 '333']
particle_name_list = ['charged_hadron', 'pion_p', 'kaon_p', 'proton',
                      'pion_m', 'kaon_m', 'anti_proton',
                      'Lambda', 'anti_Lambda', 'Xi_m', 'anti_Xi_p',
                      'Omega', 'anti_Omega', 'phi']

n_order = 7


def calculate_meanpT_inte_vn(pT_low, pT_high, data):
    """
        this function calculates the dN/dy, <pT>, and pT-integrated vn in a
        given pT range (pT_low, pT_high) for every event in the data
    """
    dN_event = data[:, 2]
    pT_event = data[:, 0]

    pT_inte_array = linspace(0., 4.0, 200)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    N = sum(dN_interp*pT_inte_array)*dpT*2.*pi

    npT = 50
    pT_inte_array = linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    meanpT = sum(dN_interp*pT_inte_array**2.)/sum(dN_interp*pT_inte_array)
    temp_vn_array = [N, meanpT,]
    for iorder in range(1, n_order):
        vn_real_event = data[:, 4*iorder]
        vn_imag_event = data[:, 4*iorder+2]
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


def calculate_rhon(data_arr):
    """
        this function calculates the rho_n correlator between vn and <pT>
        rho_n = (<\hat{delta} vn^2  \hat{delta} pT>)
                /sqrt{<(\hat{delta} vn^2)^2><\hat{delta} pT^2>}
        \hat{delta}O = delta O - delta O*delta N/<(N - <N>)^2>*delta N

        data_arr = [N, <pT>, vn]
    """
    nev = len(data_arr[:, 0])
    meanN = real(mean(data_arr[:, 0]))
    varN  = std(data_arr[:, 0])**2.

    delta_N  = data_arr[:, 0] - meanN
    delta_pT = data_arr[:, 1] - mean(data_arr[:, 1])
    delta_v2 = abs(data_arr[:, 3])**2. - mean(abs(data_arr[:, 3])**2.)
    delta_v3 = abs(data_arr[:, 4])**2. - mean(abs(data_arr[:, 4])**2.)
    delta_v4 = abs(data_arr[:, 5])**2. - mean(abs(data_arr[:, 5])**2.)

    hat_delta_pT = delta_pT - mean(delta_pT*delta_N)/varN*delta_N
    hat_delta_v2 = delta_v2 - mean(delta_v2*delta_N)/varN*delta_N
    hat_delta_v3 = delta_v3 - mean(delta_v3*delta_N)/varN*delta_N
    hat_delta_v4 = delta_v4 - mean(delta_v4*delta_N)/varN*delta_N

    # compute the error using jack-knife
    rho2_array = zeros(nev)
    rho3_array = zeros(nev)
    rho4_array = zeros(nev)
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        rho2_array[iev] = (mean(hat_delta_v2[array_idx]*hat_delta_pT[array_idx])
                           /sqrt(mean(hat_delta_v2[array_idx]**2.)
                                 *mean(hat_delta_pT**2.)))
        rho3_array[iev] = (mean(hat_delta_v3[array_idx]*hat_delta_pT[array_idx])
                           /sqrt(mean(hat_delta_v3[array_idx]**2.)
                                 *mean(hat_delta_pT**2.)))
        rho4_array[iev] = (mean(hat_delta_v4[array_idx]*hat_delta_pT[array_idx])
                           /sqrt(mean(hat_delta_v4[array_idx]**2.)
                                 *mean(hat_delta_pT**2.)))
    rho2_mean  = real(mean(rho2_array))
    rho2_err   = real(sqrt((nev - 1.)/nev*sum((rho2_array - rho2_mean)**2.)))
    rho3_mean  = real(mean(rho3_array))
    rho3_err   = real(sqrt((nev - 1.)/nev*sum((rho3_array - rho3_mean)**2.)))
    rho4_mean  = real(mean(rho4_array))
    rho4_err   = real(sqrt((nev - 1.)/nev*sum((rho4_array - rho4_mean)**2.)))
    return(array([meanN, rho2_mean, rho2_err, rho3_mean, rho3_err,
                  rho4_mean, rho4_err]))

def calculate_rhonm(data_arr):
    """
        this function calculates the rho_n correlator between vn and <pT>
        rho_n = (<\hat{delta} vn^2  \hat{delta} pT>)
                /sqrt{<(\hat{delta} vn^2)^2><\hat{delta} pT^2>}
        \hat{delta}O = delta O - delta O*delta N/<(N - <N>)^2>*delta N

        data_arr = [N, <pT>, vn]
    """
    nev = len(data_arr[:, 0])
    meanN = mean(data_arr[:, 0])
    varN  = std(data_arr[:, 0])**2.

    delta_N  = data_arr[:, 0] - meanN
    delta_pT = data_arr[:, 1] - mean(data_arr[:, 1])
    delta_v2 = abs(data_arr[:, 3])**2. - mean(abs(data_arr[:, 3])**2.)
    delta_v3 = abs(data_arr[:, 4])**2. - mean(abs(data_arr[:, 4])**2.)
    delta_v4 = abs(data_arr[:, 5])**2. - mean(abs(data_arr[:, 5])**2.)

    hat_delta_pT = delta_pT - mean(delta_pT*delta_N)/varN*delta_N
    hat_delta_v2 = delta_v2 - mean(delta_v2*delta_N)/varN*delta_N
    hat_delta_v3 = delta_v3 - mean(delta_v3*delta_N)/varN*delta_N
    hat_delta_v4 = delta_v4 - mean(delta_v4*delta_N)/varN*delta_N

    # compute the error using jack-knife
    rho23_array = zeros(nev)
    rho24_array = zeros(nev)
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        rho23_array[iev] = (
            mean(hat_delta_v2[array_idx]*hat_delta_v3[array_idx]
                 *hat_delta_pT[array_idx])
            /sqrt(mean(hat_delta_v2[array_idx]**2.)
                  *mean(hat_delta_v3[array_idx]**2.)
                  *mean(hat_delta_pT**2.)))
        rho24_array[iev] = (
            mean(hat_delta_v2[array_idx]*hat_delta_v4[array_idx]
                 *hat_delta_pT[array_idx])
            /sqrt(mean(hat_delta_v2[array_idx]**2.)
                  *mean(hat_delta_v4[array_idx]**2.)
                  *mean(hat_delta_pT**2.)))
    rho23_mean  = mean(rho23_array)
    rho23_err   = sqrt((nev - 1.)/nev*sum((rho23_array - rho23_mean)**2.))
    rho24_mean  = mean(rho24_array)
    rho24_err   = sqrt((nev - 1.)/nev*sum((rho24_array - rho24_mean)**2.))
    return(array([meanN, rho23_mean, rho23_err, rho24_mean, rho24_err]))


def calculate_meanpT_moments(data_arr):
    """
        data_arr = [N, <pT>, vn]
    """
    nev = len(data_arr[:, 0])
    meanN = real(mean(data_arr[:, 0]))
    meanpT = real(mean(data_arr[:, 1]))
    meanpT_err = real(std(data_arr[:, 1]))/sqrt(nev)
    varN  = std(data_arr[:, 0])**2.

    delta_N  = data_arr[:, 0] - meanN
    delta_pT = data_arr[:, 1] - mean(data_arr[:, 1])

    hat_delta_pT = delta_pT - mean(delta_pT*delta_N)/varN*delta_N

    # compute the error using jack-knife
    rho2_array = zeros(nev)
    rho3_array = zeros(nev)
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        rho2_array[iev] = (sqrt(mean(hat_delta_pT[array_idx]**2.))
                           /mean(data_arr[array_idx, 1]))
        rho3_array[iev] = (mean(hat_delta_pT[array_idx]**3.)
                           /(mean(hat_delta_pT[array_idx]**2.)**1.5))
    rho2_mean  = real(mean(rho2_array))
    rho2_err   = real(sqrt((nev - 1.)/nev*sum((rho2_array - rho2_mean)**2.)))
    rho3_mean  = real(mean(rho3_array))
    rho3_err   = real(sqrt((nev - 1.)/nev*sum((rho3_array - rho3_mean)**2.)))
    return(array([meanN, meanpT, meanpT_err, rho2_mean, rho2_err,
                  rho3_mean, rho3_err]))


def calculate_meanpT_fluc(dN_array, pT_array, pT_min=0.0, pT_max=3.0):
    """
        This function computes the mean pT fluctuations
            returns sigma_pT/<pT>, sigma_pT/<pT>_err
        here sigma_pT is the standard deviation of the event-by-event mean pT
        This function accepts pT_cut through [pT_min, pT_max]
        dN_array is dN/(2\pi dy pT dpT)
    """
    npT_interp = 50
    pT_inte_array = linspace(pT_min, pT_max, npT_interp)

    nev, npT = dN_array.shape
    mean_pT_array = zeros(nev)
    for iev in range(nev):
        dN_interp = exp(interp(pT_inte_array, pT_array[iev, :],
                               log(dN_array[iev, :] + 1e-30)))
        mean_pT_array[iev] = (sum(pT_inte_array**2.*dN_interp)
                              /sum(pT_inte_array*dN_interp))

    # compute the error using jack-knife
    rn_array  = zeros(nev)
    for iev in range(nev):
        array_idx      = [True]*nev
        array_idx[iev] = False
        array_idx      = array(array_idx)
        rn_ev          = (std(mean_pT_array[array_idx])
                          /(mean(mean_pT_array[array_idx]) + 1e-15))
        rn_array[iev]  = rn_ev
    rn_mean  = mean(rn_array, axis=0)
    rn_err   = sqrt((nev - 1.)/nev*sum((rn_array - rn_mean)**2.))
    return([rn_mean, rn_err])


def calculate_symmetric_cumulant(data_arr):
    """
        this funciton computes the symmetric cumulant
            SC(m,n) = <v_m*conj(v_m)*v_n*conj(v_n)>
                      - <v_m*conj(v_m)>*<v_n*conj(v_n)>
            NSC(m,n) = (<v_m*conj(v_m)*v_n*conj(v_n)>
                        /<v_m*conj(v_m)>*<v_n*conj(v_n)> - 1)
        we use Jackknife resampling method to estimate the statistical error
        data_arr = [N, <pT>, vn]
        return(SC, NSC) for (2,3) and (2,4)
    """
    nev = len(data_arr[:, 0])
    dN = real(data_arr[:, 0])
    meanN = mean(dN)
    dN = real(data_arr[:, -1])
    Q1 = dN*data_arr[:, 2]
    Q2 = dN*data_arr[:, 3]
    Q3 = dN*data_arr[:, 4]
    Q4 = dN*data_arr[:, 5]
    Q5 = dN*data_arr[:, 6]
    Q6 = dN*data_arr[:, 7]

    # two-particle correlation
    N2_weight = dN*(dN - 1.)
    Q2_2 = abs(Q2)**2. - dN
    Q3_2 = abs(Q3)**2. - dN
    Q4_2 = abs(Q4)**2. - dN

    # four-particle correlation
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q_32 = ((abs(Q2)**2.)*(abs(Q3)**2.) - 2.*real(Q5*conj(Q2)*conj(Q3))
        - 2.*real(Q3*conj(Q1)*conj(Q2)) + abs(Q5)**2. + abs(Q1)**2.
        - (dN - 4.)*(abs(Q2)**2. + abs(Q3)**2.) + dN*(dN - 6.)
    )
    Q_42 = ((abs(Q2)**2.)*(abs(Q4)**2.) - 2.*real(Q6*conj(Q2)*conj(Q4))
        - 2.*real(Q4*conj(Q2)*conj(Q2)) + abs(Q6)**2. + abs(Q2)**2.
        - (dN - 4.)*(abs(Q2)**2. + abs(Q4)**2.) + dN*(dN - 6.)
    )

    # calcualte observables with Jackknife resampling method
    SC32_array = zeros(nev)
    SC42_array = zeros(nev)
    NSC32_array = zeros(nev)
    NSC42_array = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        # SC(3,2)
        v2v3 = ((mean(Q3_2[array_idx])*mean(Q2_2[array_idx]))
                /(mean(N2_weight[array_idx])**2.))

        SC32_array[iev] = (
                mean(Q_32[array_idx])/mean(N4_weight[array_idx]) - v2v3)
        NSC32_array[iev] = SC32_array[iev]/v2v3

        # SC(4,2)
        v2v4 = ((mean(Q4_2[array_idx])*mean(Q2_2[array_idx]))
                /(mean(N2_weight[array_idx])**2.))
        SC42_array[iev] = (
                mean(Q_42[array_idx])/mean(N4_weight[array_idx]) - v2v4)
        NSC42_array[iev] = SC42_array[iev]/v2v4

    SC32_mean = mean(SC32_array)
    SC32_err = sqrt((nev - 1.)/nev*sum((SC32_array - SC32_mean)**2.))
    NSC32_mean = mean(NSC32_array)
    NSC32_err = sqrt((nev - 1.)/nev*sum((NSC32_array - NSC32_mean)**2.))

    SC42_mean = mean(SC42_array)
    SC42_err = sqrt((nev - 1.)/nev*sum((SC42_array - SC42_mean)**2.))
    NSC42_mean = mean(NSC42_array)
    NSC42_err = sqrt((nev - 1.)/nev*sum((NSC42_array - NSC42_mean)**2.))

    results = [meanN, SC32_mean, SC32_err, SC42_mean, SC42_err,
               NSC32_mean, NSC32_err, NSC42_mean, NSC42_err]
    return(results)


hf = h5py.File(data_path, "r")
event_list = list(hf.keys())
print("total number of events: {}".format(len(event_list)))

dN_dy_mb = zeros(len(event_list))
for ifolder, event_name in enumerate(event_list):
    file_name   = "particle_9999_vndata_eta_-0.5_0.5.dat"
    try:
        event_group = hf.get(event_name)
        temp_data   = event_group.get(file_name)
        temp_data   = nan_to_num(temp_data)
        dN_dy_mb[ifolder] = -temp_data[0, 1]
    except:
        continue
dN_dy_mb = -sort(dN_dy_mb)


output_filename = "charged_hadron_meanpT_fluct_ALICE.dat"
ch_pT_file_alice = open(path.join(avg_folder_header, output_filename), 'w')
ch_pT_file_alice.write(
        "# dN/deta  <pT>  <pT>_err  <(delta pT)^2>  <(delta pT)^2>_err  "
        "<(delta pT)^3>   <(delta pT)^3>_err\n")
output_filename = "charged_hadron_symmetric_cumulants_ALICE.dat"
ch_SC_file_alice = open(path.join(avg_folder_header, output_filename), 'w')
ch_SC_file_alice.write(
        "# dN/deta  SC{2,3}  SC{2,3}_err  SC{2,4}  SC{2,4}_err "
        + "NSC{2,3}  NSC{2,3}_err  NSC{2,4}  NSC{2,4}_err\n")
output_filename = "charged_hadron_rho_n_ALICE.dat"
ch_rho_file_alice = open(path.join(avg_folder_header, output_filename), 'w')
ch_rho_file_alice.write(
        "# dN/deta  rho_2  rho_2_err  rho_3  rho_3_err  rho_4  rho_4_err\n")
output_filename = "charged_hadron_rho_nm_ALICE.dat"
ch_rho_nm_file_alice = open(path.join(avg_folder_header, output_filename), 'w')
ch_rho_nm_file_alice.write(
        "# dN/deta  rho_23  rho_23_err  rho_24  rho_24_err\n")
output_filename = "charged_hadron_rho_n_ATLAS.dat"
ch_rho_file_atlas = open(path.join(avg_folder_header, output_filename), 'w')
ch_rho_file_atlas.write(
        "# dN/deta  rho_2  rho_2_err  rho_3  rho_3_err  rho_4  rho_4_err\n")
output_filename = "charged_hadron_rho_nm_ATLAS.dat"
ch_rho_nm_file_atlas = open(path.join(avg_folder_header, output_filename), 'w')
ch_rho_nm_file_atlas.write(
        "# dN/deta  rho_23  rho_23_err  rho_24  rho_24_err\n")
for icen in range(len(centrality_cut_list) - 1):
    if centrality_cut_list[icen+1] < centrality_cut_list[icen]: continue
    avg_folder = path.join(
        avg_folder_header, "{0:02.0f}-{1:02.0f}".format(
            centrality_cut_list[icen], centrality_cut_list[icen+1])
    )
    mkdir(avg_folder)

    selected_events_list = []
    for ifolder, event_name in enumerate(event_list):
        dN_dy_cut_high = (
            dN_dy_mb[int(len(dN_dy_mb)*centrality_cut_list[icen]/100.)])
        dN_dy_cut_low  = dN_dy_mb[
            min(len(dN_dy_mb)-1,
                int(len(dN_dy_mb)*centrality_cut_list[icen+1]/100.))
        ]
        file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
        try:
            event_group = hf.get(event_name)
            temp_data   = event_group.get(file_name)
            temp_data   = nan_to_num(temp_data)
            if (temp_data[0, 1] > dN_dy_cut_low
                and temp_data[0, 1] <= dN_dy_cut_high):
                selected_events_list.append(event_name)
        except:
            continue

    nev = len(selected_events_list)
    print("analysis {}%-{}% nev = {}...".format(
            centrality_cut_list[icen], centrality_cut_list[icen+1], nev))
    if nev == 0:
        print("Skip ...")
        continue

    for ipart, particle_id in enumerate([particle_list[0]]):
        print("processing %s ..." % particle_name_list[ipart])

        # first particle yield dN/dy
        if particle_id == '9999':
            file_name = 'particle_%s_vndata_eta_-0.5_0.5.dat' % particle_id
        else:
            file_name = 'particle_%s_vndata_y_-0.5_0.5.dat' % particle_id

        dN_dy = []
        for ifolder, event_name in enumerate(selected_events_list):
            event_group = hf.get(event_name)
            temp_data = event_group.get(file_name)
            temp_data = nan_to_num(temp_data)

            dN_dy.append(temp_data[0, 1])

        dN_dy = array(dN_dy)
        dN_dy_avg = mean(dN_dy)
        dN_dy_avg_err = std(dN_dy)/sqrt(nev)

        # then <pT>, vn, dN/(2pi dy pT dpT), vn{SP}(pT)
        if particle_id == '9999':
            #file_name = 'particle_%s_vndata_diff_eta_-0.5_0.5.dat' % particle_id
            file_name = 'particle_%s_vndata_diff_eta_-0.8_0.8.dat' % particle_id
        else:
            file_name = 'particle_%s_vndata_diff_y_-0.5_0.5.dat' % particle_id
        file_name_ref1 = 'particle_9999_vndata_diff_eta_0.5_2.dat'
        file_name_ref2 = 'particle_9999_vndata_diff_eta_-2_-0.5.dat'
        file_name_ALICE = 'particle_9999_vndata_diff_eta_-0.8_0.8.dat'
        file_name_ATLAS = 'particle_9999_vndata_diff_eta_-2_2.dat'
        ecc_filename = "eccentricities_evo_eta_-0.5_0.5.dat"
        eccp_filename = "momentum_anisotropy_eta_-0.5_0.5.dat"

        pT_array = []
        dN_array = []
        pT_array_ALICE = []
        dN_array_ALICE = []
        pT_array_ATLAS = []
        dN_array_ATLAS = []
        vn_phenix_array = []
        vn_phenix_array_sub1 = []; vn_phenix_array_sub2 = []
        vn_star_array = []
        vn_star_array_sub1 = []; vn_star_array_sub2 = []
        vn_alice_array = []
        vn_alice_array_sub1 = []; vn_alice_array_sub2 = []
        vn_cms_array = []
        vn_cms_array_sub1 = []; vn_cms_array_sub2 = []
        vn_cms_arrays_for_rn = []
        vn_atlas_array = []
        vn_atlas_array_sub1 = []; vn_atlas_array_sub2 = []
        eccn_array = []
        eccp2_array = []
        for ifolder, event_name in enumerate(selected_events_list):
            event_group = hf.get(event_name)
            temp_data = event_group.get(file_name)
            temp_data = nan_to_num(temp_data)
            temp_data_ref1 = event_group.get(file_name_ref1)
            temp_data_ref1 = nan_to_num(temp_data_ref1)
            temp_data_ref2 = event_group.get(file_name_ref2)
            temp_data_ref2 = nan_to_num(temp_data_ref2)
            temp_ecc_data = event_group.get(ecc_filename)
            temp_eccp_data = event_group.get(eccp_filename)

            dN_event = temp_data[:, 2]  # dN/(2pi dy pT dpT)
            pT_event = temp_data[:, 0]

            init_eccn = [
                temp_ecc_data[0, 2*i-1] + 1j*temp_ecc_data[0, 2*i] for i in range(1, n_order)]
            init_eccp = [
                temp_eccp_data[0, 2*i-1] + 1j*temp_eccp_data[0, 2*i] for i in range(1, 4)]

            # record particle spectra
            pT_array.append(pT_event)
            dN_array.append(dN_event)
            eccn_array.append(init_eccn)
            eccp2_array.append(init_eccp)
            if particle_id == "9999":
                try:
                    temp_data1 = event_group.get(file_name_ATLAS)
                    temp_data1 = nan_to_num(temp_data1)
                    pT_array_ATLAS.append(temp_data1[:, 0])
                    dN_array_ATLAS.append(temp_data1[:, 2])
                    temp_data1 = event_group.get(file_name_ALICE)
                    temp_data1 = nan_to_num(temp_data1)
                    pT_array_ALICE.append(temp_data1[:, 0])
                    dN_array_ALICE.append(temp_data1[:, 2])
                except:
                    continue

            # pT-integrated vn
            # vn with PHENIX pT cut
            temp_vn_array = calculate_meanpT_inte_vn(0.2, 2.0, temp_data)
            vn_phenix_array.append(temp_vn_array)

            # vn with STAR pT cut
            temp_vn_array = calculate_meanpT_inte_vn(0.15, 2.0, temp_data)
            vn_star_array.append(temp_vn_array)

            # vn with ALICE pT cut
            temp_vn_array = calculate_meanpT_inte_vn(0.2, 3.0, temp_data)
            vn_alice_array.append(temp_vn_array)

            # vn with CMS pT cut
            temp_vn_array = calculate_meanpT_inte_vn(0.3, 3.0, temp_data)
            vn_cms_array.append(temp_vn_array)

            # vn with ATLAS pT cut
            temp_vn_array = calculate_meanpT_inte_vn(0.5, 3.0, temp_data)
            vn_atlas_array.append(temp_vn_array)

        # now we perform event average
        dN_array = array(dN_array)
        pT_array = array(pT_array)

        vn_phenix_array = array(vn_phenix_array)
        vn_star_array   = array(vn_star_array)
        vn_alice_array  = array(vn_alice_array)
        vn_cms_array    = array(vn_cms_array)
        vn_atlas_array  = array(vn_atlas_array)

        #rhon_phenix = calculate_rhon(vn_phenix_array)
        #rhon_star   = calculate_rhon(vn_star_array)
        #rhon_cms    = calculate_rhon(vn_cms_array)
        rhon_atlas   = calculate_rhon(vn_atlas_array)
        rhonm_atlas  = calculate_rhonm(vn_atlas_array)
        rhon_alice  = calculate_rhon(vn_alice_array)
        rhonm_alice  = calculate_rhonm(vn_alice_array)
        SC_alice = calculate_symmetric_cumulant(vn_alice_array)
        meanpT_Cn_alice = calculate_meanpT_moments(vn_alice_array)
        #print(meanpT_Cn_alice)

        n_pT = len(pT_array[0, :])
        pT_spectra = zeros([n_pT])
        for ipT in range(len(pT_array[0, :])):
            dN_temp = sum(dN_array[:, ipT]*pT_array[:, ipT])
            if(dN_temp > 0):
                pT_spectra[ipT] = (
                        sum(pT_array[:, ipT]**2.*dN_array[:, ipT])/dN_temp)
            else:
                pT_spectra[ipT] = mean(pT_array[:, ipT])
        # dN/(2pi dy pT dpT)
        dN_spectra = mean(pT_array*dN_array, 0)/pT_spectra
        dN_spectra_err = std(pT_array*dN_array, 0)/pT_spectra/sqrt(nev)

        # compute mean pT event-by-event
        sigma_pT, sigma_pT_err = calculate_meanpT_fluc(dN_array, pT_array,
                                                       0.15, 2.0)

        if particle_id == "9999":
            # calculate dNch (pT > 0.4 |eta| < 2.5) for ATLAS
            dN_array_ATLAS = array(dN_array_ATLAS)
            pT_array_ATLAS = array(pT_array_ATLAS)

            n_pT = len(pT_array_ATLAS[0, :])
            pT_spectra_ATLAS = zeros([n_pT])
            for ipT in range(len(pT_array_ATLAS[0, :])):
                dN_temp = sum(dN_array_ATLAS[:, ipT]*pT_array_ATLAS[:, ipT])
                if (dN_temp > 0):
                    pT_spectra_ATLAS[ipT] = (
                        sum(pT_array_ATLAS[:, ipT]**2.*dN_array_ATLAS[:, ipT])
                        /dN_temp)
                else:
                    pT_spectra_ATLAS[ipT] = mean(pT_array_ATLAS[:, ipT])
            # dN/(2pi pT dpT) (ATLAS)
            dN_spectra_ATLAS = (mean(pT_array_ATLAS*dN_array_ATLAS, 0)
                                /(pT_spectra_ATLAS + 1e-30))
            dN_spectra_ATLAS_err = (std(pT_array_ATLAS*dN_array_ATLAS, 0)
                                    /(pT_spectra_ATLAS + 1e-30)/sqrt(nev))
            pT_interp1 = linspace(0.4, 3.0, 30)
            dpT1 = pT_interp1[1] - pT_interp1[0]
            dN_interp1 = exp(interp(pT_interp1, pT_spectra_ATLAS,
                                    log(dN_spectra_ATLAS + 1e-30)))
            dN_interp1_err = interp(pT_interp1, pT_spectra_ATLAS,
                                    dN_spectra_ATLAS_err)
            dNch_ATLAS = sum(pT_interp1*dN_interp1)*dpT1*2*pi*5.
            print("dNch(pT > 0.4 GeV, |eta| < 2.5) = {0:.4f}".format(
                                                                dNch_ATLAS))

            # calculate dNch (pT > 0.2 |eta| < 0.8) for ALICE
            dN_array_ALICE = array(dN_array_ALICE)
            pT_array_ALICE = array(pT_array_ALICE)

            n_pT = len(pT_array_ALICE[0, :])
            pT_spectra_ALICE = zeros([n_pT])
            for ipT in range(len(pT_array_ALICE[0, :])):
                dN_temp = sum(dN_array_ALICE[:, ipT]*pT_array_ALICE[:, ipT])
                if (dN_temp > 0):
                    pT_spectra_ALICE[ipT] = (
                        sum(pT_array_ALICE[:, ipT]**2.*dN_array_ALICE[:, ipT])
                        /dN_temp)
                else:
                    pT_spectra_ALICE[ipT] = mean(pT_array_ALICE[:, ipT])
            # dN/(2pi pT dpT)
            dN_spectra_ALICE = (mean(pT_array_ALICE*dN_array_ALICE, 0)
                                /(pT_spectra_ALICE + 1e-30))
            dN_spectra_ALICE_err = (std(pT_array_ALICE*dN_array_ALICE, 0)
                                    /(pT_spectra_ALICE + 1e-30)/sqrt(nev))
            pT_interp1 = linspace(0.2, 3.0, 30)
            dpT1 = pT_interp1[1] - pT_interp1[0]
            dN_interp1 = exp(interp(pT_interp1, pT_spectra_ALICE,
                                    log(dN_spectra_ALICE + 1e-30)))
            dN_interp1_err = interp(pT_interp1, pT_spectra_ALICE,
                                    dN_spectra_ALICE_err)
            dNch_ALICE = sum(pT_interp1*dN_interp1)*dpT1*2*pi*1.6
            print("dNch(pT > 0.2 GeV, |eta| < 0.8) = {0:.4f}".format(
                                                                dNch_ALICE))

            # calculate dNch (pT > 0.15 |eta| < 0.8) for ALICE
            pT_interp1 = linspace(0.15, 3.0, 30)
            dpT1 = pT_interp1[1] - pT_interp1[0]
            dN_interp1 = exp(interp(pT_interp1, pT_spectra_ALICE,
                                    log(dN_spectra_ALICE + 1e-30)))
            dN_interp1_err = interp(pT_interp1, pT_spectra_ALICE,
                                    dN_spectra_ALICE_err)
            dNch_ALICE2 = sum(pT_interp1*dN_interp1)*dpT1*2*pi*1.6
            print("dNch(pT > 0.15 GeV, |eta| < 0.8) = {0:.4f}".format(
                                                                dNch_ALICE2))

        # calculate mean pT from event-averaged particle spectrum
        pT_interp = linspace(0.05, 2.95, 30)
        dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra + 1e-30)))
        dN_interp_err = interp(pT_interp, pT_spectra, dN_spectra_err)
        mean_pT = sum(pT_interp**2.*dN_interp)/sum(pT_interp*dN_interp)
        mean_pT_upper = (sum(pT_interp**2.*(dN_interp+dN_interp_err))
                         /sum(pT_interp*(dN_interp+dN_interp_err)))
        mean_pT_lower = (sum(pT_interp**2.*(dN_interp-dN_interp_err))
                         /sum(pT_interp*(dN_interp-dN_interp_err)))
        mean_pT_err = max(abs(mean_pT_upper - mean_pT), 
                          abs(mean_pT - mean_pT_lower))
        pT_interp = linspace(0.15, 2.95, 30)
        dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra + 1e-30)))
        dN_interp_err = interp(pT_interp, pT_spectra, dN_spectra_err)
        mean_pT_1 = sum(pT_interp**2.*dN_interp)/sum(pT_interp*dN_interp)
        mean_pT_1_upper = (sum(pT_interp**2.*(dN_interp+dN_interp_err))
                         /sum(pT_interp*(dN_interp+dN_interp_err)))
        mean_pT_1_lower = (sum(pT_interp**2.*(dN_interp-dN_interp_err))
                         /sum(pT_interp*(dN_interp-dN_interp_err)))
        mean_pT_1_err = max(abs(mean_pT_1_upper - mean_pT_1), 
                          abs(mean_pT_1 - mean_pT_1_lower))

        ######################################################################
        # finally, output all the results
        ######################################################################
        vn_alice_array = array(vn_alice_array)
        if particle_id == "9999":
            eccn_array = array(eccn_array)
            eccp2_array = array(eccp2_array)
            for iorder in range(2, 4):
                output_filename = ("ebe_v{}_ALICE.dat".format(iorder))
                f = open(path.join(avg_folder, output_filename), 'w')
                if iorder == 2:
                    f.write("# ecc_n(real) ecc_n(imag) v_n(real) v_n(imag) "
                            + "ecc_p2(real) ecc_p2(imag)\n")
                else:
                    f.write("# ecc_n(real) ecc_n(imag) v_n(real) v_n(imag)\n")
                for iev in range(len(eccn_array[:, 0])):
                    if iorder == 2:
                        f.write("%.5e  %.5e  %.5e  %.5e  %.5e  %.5e\n"
                                % (real(eccn_array[iev, iorder-1]),
                                   imag(eccn_array[iev, iorder-1]),
                                   real(vn_alice_array[iev, iorder+1]),
                                   imag(vn_alice_array[iev, iorder+1]),
                                   real(eccp2_array[iev, 2]),
                                   imag(eccp2_array[iev, 2]))
                        )
                    else:
                        f.write("%.5e  %.5e  %.5e  %.5e\n"
                                % (real(eccn_array[iev, iorder-1]),
                                   imag(eccn_array[iev, iorder-1]),
                                   real(vn_alice_array[iev, iorder+1]),
                                   imag(vn_alice_array[iev, iorder+1]))
                        )
                f.close()

        # output rho_n
        if particle_id == "9999":
            for val in rhon_alice:
                ch_rho_file_alice.write("%.5e  " % val)
            ch_rho_file_alice.write("\n")
            for val in rhonm_alice:
                ch_rho_nm_file_alice.write("%.5e  " % val)
            ch_rho_nm_file_alice.write("\n")
            for val in rhon_atlas:
                ch_rho_file_atlas.write("%.5e  " % val)
            ch_rho_file_atlas.write("\n")
            for val in rhonm_atlas:
                ch_rho_nm_file_atlas.write("%.5e  " % val)
            ch_rho_nm_file_atlas.write("\n")
            for val in SC_alice:
                ch_SC_file_alice.write("%.5e  " % val)
            ch_SC_file_alice.write("\n")
            for val in meanpT_Cn_alice:
                ch_pT_file_alice.write("%.5e  " % val)
            ch_pT_file_alice.write("\n")

        output_filename = ("%s_integrated_observables.dat"
                           % particle_name_list[ipart])
        f = open(path.join(avg_folder, output_filename), 'w')
        f.write("dN/dy= %.5e +/- %.5e\n" % (dN_dy_avg, dN_dy_avg_err))
        if particle_id == "9999":
            f.write("dN/dy(pT>0.2,|eta|<0.8)= %.5e \n" % (dNch_ALICE))
            f.write("dN/dy(pT>0.15,|eta|<0.8)= %.5e \n" % (dNch_ALICE2))
            f.write("dN/dy(pT>0.4,|eta|<2.5)= %.5e \n" % (dNch_ATLAS))
        f.write("<pT>= %.5e +/- %.5e\n" % (mean_pT, mean_pT_err))
        f.write("<pT(>0.15)>= %.5e +/- %.5e\n" % (mean_pT_1, mean_pT_1_err))
        f.write("sigma_pT/<pT>= %.5e +/- %.5e\n"% (sigma_pT, sigma_pT_err))
        f.close()

        output_filename = (
            "{}_ebe_N_meanpT_vn_ALICE.dat".format(particle_name_list[ipart]))
        f = open(path.join(avg_folder, output_filename), 'w')
        f.write("# dN/dy  <pT>(GeV)  v_n(real)  v_n(imag)\n")
        for iev in range(len(vn_alice_array[:, 0])):
            f.write("%.5e  %.5e  " % (real(vn_alice_array[iev, 0]),
                                      real(vn_alice_array[iev, 1])))
            for iorder in range(1, n_order):
                f.write("%.5e  %.5e  " % (real(vn_alice_array[iev, iorder+1]),
                                          imag(vn_alice_array[iev, iorder+1])))
            f.write("\n")
        f.close()

        output_filename = (
            "{}_ebe_N_meanpT_vn_ATLAS.dat".format(particle_name_list[ipart]))
        f = open(path.join(avg_folder, output_filename), 'w')
        f.write("# dN/dy  <pT>(GeV)  v_n(real)  v_n(imag)\n")
        for iev in range(len(vn_atlas_array[:, 0])):
            f.write("%.5e  %.5e  " % (real(vn_atlas_array[iev, 0]),
                                      real(vn_atlas_array[iev, 1])))
            for iorder in range(1, n_order):
                f.write("%.5e  %.5e  " % (real(vn_atlas_array[iev, iorder+1]),
                                          imag(vn_atlas_array[iev, iorder+1])))
            f.write("\n")
        f.close()

ch_rho_file_atlas.close()
ch_rho_nm_file_atlas.close()
ch_rho_file_alice.close()
ch_rho_nm_file_alice.close()
ch_SC_file_alice.close()
ch_pT_file_alice.close()
print("Analysis is done.")

