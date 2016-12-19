#! /usr/bin/env python2.7
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
import shutil

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

try:
    working_folder = path.abspath(argv[1])
    avg_folder = path.join(path.abspath(argv[2]),
                           working_folder.split('/')[-1])
    print("output folder: %s" % avg_folder)
    if(path.isdir(avg_folder)):
        print("folder %s already exists!" % avg_folder)
        var = raw_input("do you want to delete it? [y/N]")
        if 'y' in var:
            shutil.rmtree(avg_folder)
        else:
            print("please choose another folder path~")
            exit(0)
    mkdir(avg_folder)
except IndexError:
    print("Usage: average_event_spvn.py working_folder results_folder")
    exit(1)

particle_list = ['pion_p', 'pion_m', 'Kaon_p', 'Kaon_m', 'proton',
                 'anti_proton', 'Lambda', 'Xi_m', 'anti_Xi',
                 'Omega', 'Charged_eta']
nonlinear_reponse_correlator_name_list = [
                'v4_L', 'v4(Psi2)', 'rho_422', 'chi_422',
                'v5_L', 'v5(Psi23)', 'rho_523', 'chi_523',
                'v6_L', 'v6(Psi2)', 'v6(Psi3)',
                'rho_6222', 'rho_633', 'chi_6222', 'chi_633']

n_order = 7

def calcualte_inte_vn(pT_low, pT_high, data):
    npT = 50
    pT_inte_array = linspace(pT_low, pT_high, npT)
    dN_event = data[:, 2]
    pT_event = data[:, 0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    temp_vn_array = []
    for iorder in range(1, n_order):
        vn_real_event = data[:, 3*iorder]
        vn_imag_event = data[:, 3*iorder+1]
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


def calculate_chi_422(vn_array):
    v2_array = vn_array[:, 1]
    nev = len(v2_array)
    v4_array = vn_array[:, 3]
    chi_422_num = v4_array*(v2_array.conjugate())**2
    chi_422_den = (abs(v2_array))**4
    chi_422_num_ave = mean(chi_422_num)
    chi_422_den_ave = mean(chi_422_den)
    chi_422_num_err = std(chi_422_num)/sqrt(nev)
    chi_422_den_err = std(chi_422_den)/sqrt(nev)
    chi_422 = chi_422_num_ave/chi_422_den_ave
    chi_422_err = sqrt(
        (chi_422_num_err/chi_422_den_ave)**2.
        + (chi_422_num_ave*chi_422_den_err/(chi_422_den_ave)**2.)**2.)
    return(chi_422.real, chi_422_err.real)


def calculate_chi_523(vn_array):
    v2_array = vn_array[:, 1]
    nev = len(v2_array)
    v3_array = vn_array[:, 2]
    v5_array = vn_array[:, 4]
    chi_523_num = v5_array*(v2_array.conjugate()*v3_array.conjugate())
    chi_523_den = (abs(v2_array))**2*(abs(v3_array))**2
    chi_523_num_ave = mean(chi_523_num)
    chi_523_den_ave = mean(chi_523_den)
    chi_523_num_err = std(chi_523_num)/sqrt(nev)
    chi_523_den_err = std(chi_523_den)/sqrt(nev)
    chi_523 = chi_523_num_ave/chi_523_den_ave
    chi_523_err = sqrt(
        (chi_523_num_err/chi_523_den_ave)**2.
        + (chi_523_num_ave*chi_523_den_err/(chi_523_den_ave)**2.)**2.)
    return(chi_523.real, chi_523_err.real)


def calculate_chi_6222(vn_array):
    v6_array = vn_array[:, 5]
    v2_array = vn_array[:, 1]
    nev = len(v2_array)
    chi_6222_num = v6_array*((v2_array.conjugate())**3)
    chi_6222_den = (abs(v2_array))**6
    chi_6222_num_ave = mean(chi_6222_num)
    chi_6222_den_ave = mean(chi_6222_den)
    chi_6222_num_err = std(chi_6222_num)/sqrt(nev)
    chi_6222_den_err = std(chi_6222_den)/sqrt(nev)
    chi_6222 = chi_6222_num_ave/chi_6222_den_ave
    chi_6222_err = sqrt(
        (chi_6222_num_err/chi_6222_den_ave)**2.
        + (chi_6222_num_ave*chi_6222_den_err/(chi_6222_den_ave)**2.)**2.)
    return(chi_6222.real, chi_6222_err.real)


def calculate_chi_633(vn_array):
    v6_array = vn_array[:, 5]
    v3_array = vn_array[:, 2]
    nev = len(v3_array)
    chi_633_num = v6_array*((v3_array.conjugate())**2)
    chi_633_den = (abs(v3_array))**4                  
    chi_633_num_ave = mean(chi_633_num)
    chi_633_den_ave = mean(chi_633_den)
    chi_633_num_err = std(chi_633_num)/sqrt(nev)
    chi_633_den_err = std(chi_633_den)/sqrt(nev)
    chi_633 = chi_633_num_ave/chi_633_den_ave
    chi_633_err = sqrt(
        (chi_633_num_err/chi_633_den_ave)**2.
        + (chi_633_num_ave*chi_633_den_err/(chi_633_den_ave)**2.)**2.)
    return(chi_633.real, chi_633_err.real)


def calculate_v4_Psi2(chi_422, chi_422_err, vn_array):
    """
        v4(Psi2) = chi_422*sqrt(<abs(V2)**4>)
    """
    v2_array = vn_array[:, 1]
    nev = len(v2_array)
    temp_array = abs(v2_array)**4
    v2_factor = sqrt(mean(temp_array))
    v2_factor_err = std(temp_array)/(2.*v2_factor)/sqrt(nev)
    v4_Psi2 = chi_422*v2_factor
    v4_Psi2_err = sqrt((chi_422_err*v2_factor)**2.
                       + (chi_422*v2_factor_err)**2.)
    return(v4_Psi2, v4_Psi2_err)


def calculate_v5_Psi23(chi_523, chi_523_err, vn_array):
    """
        v5(Psi23) = chi_523*sqrt(<abs(V2)**2*abs(V3)**2>)
    """
    v2_array = vn_array[:, 1]
    v3_array = vn_array[:, 2]
    nev = len(v2_array)
    temp_array = abs(v2_array)**2.*abs(v3_array)**2.
    v23_factor = sqrt(mean(temp_array))
    v23_factor_err = std(temp_array)/(2.*v23_factor)/sqrt(nev)
    v5_Psi23 = chi_523*v23_factor
    v5_Psi23_err = sqrt((chi_523_err*v23_factor)**2.
                        + (chi_523*v23_factor_err)**2.)
    return(v5_Psi23, v5_Psi23_err)


def calculate_v6_Psi2(chi_6222, chi_6222_err, vn_array):
    """
        v6(Psi2) = chi_6222*sqrt(<abs(V2)**6>)
    """
    v2_array = vn_array[:, 1]
    nev = len(v2_array)
    temp_array = abs(v2_array)**6.
    v2_factor = sqrt(mean(temp_array))
    v2_factor_err = std(temp_array)/(2.*v2_factor)/sqrt(nev)
    v6_Psi2 = chi_6222*v2_factor
    v6_Psi2_err = sqrt((chi_6222_err*v2_factor)**2.
                        + (chi_6222*v2_factor_err)**2.)
    return(v6_Psi2, v6_Psi2_err)


def calculate_v6_Psi3(chi_633, chi_633_err, vn_array):
    """
        v6(Psi3) = chi_633*sqrt(<abs(V3)**4>)
    """
    v3_array = vn_array[:, 2]
    nev = len(v3_array)
    temp_array = abs(v3_array)**4.
    v3_factor = sqrt(mean(temp_array))
    v3_factor_err = std(temp_array)/(2.*v3_factor)/sqrt(nev)
    v6_Psi3 = chi_633*v3_factor
    v6_Psi3_err = sqrt((chi_633_err*v3_factor)**2.
                        + (chi_633*v3_factor_err)**2.)
    return(v6_Psi3, v6_Psi3_err)


def calculate_rho_422(v4_Psi2, v4_Psi2_err, vn_array):
    """
        rho_422 = v4(Psi2)/v4(Psi4)
    """
    v4_array = vn_array[:, 3]
    nev = len(v4_array)
    v4_Psi4 = sqrt(mean(abs(v4_array)**2.))
    v4_Psi4_err = std(abs(v4_array)**2.)/(2.*v4_Psi4)/sqrt(nev)
    rho_422 = v4_Psi2/v4_Psi4
    rho_422_err = sqrt((v4_Psi2_err/v4_Psi4)**2.
                       + (v4_Psi2*v4_Psi4_err/v4_Psi4**2.)**2.)
    return(rho_422, rho_422_err)


def calculate_rho_523(v5_Psi23, v5_Psi23_err, vn_array):
    """
        rho_523 = v5(Psi23)/v5(Psi5)
    """
    v5_array = vn_array[:, 4]
    nev = len(v5_array)
    v5_Psi5 = sqrt(mean(abs(v5_array)**2.))
    v5_Psi5_err = std(abs(v5_array)**2.)/(2.*v5_Psi5)/sqrt(nev)
    rho_523 = v5_Psi23/v5_Psi5
    rho_523_err = sqrt((v5_Psi23_err/v5_Psi5)**2.
                       + (v5_Psi23*v5_Psi5_err/v5_Psi5**2.)**2.)
    return(rho_523, rho_523_err)


def calculate_rho_6222(v6_Psi2, v6_Psi2_err, vn_array):
    """
        rho_6222 = v6(Psi2)/v6(Psi6)
    """
    v6_array = vn_array[:, 5]
    nev = len(v6_array)
    v6_Psi6 = sqrt(mean(abs(v6_array)**2.))
    v6_Psi6_err = std(abs(v6_array)**2.)/(2.*v6_Psi6)/sqrt(nev)
    rho_6222 = v6_Psi2/v6_Psi6
    rho_6222_err = sqrt((v6_Psi2_err/v6_Psi6)**2.
                        + (v6_Psi2*v6_Psi6_err/v6_Psi6**2.)**2.)
    return(rho_6222, rho_6222_err)


def calculate_rho_633(v6_Psi3, v6_Psi3_err, vn_array):
    """
        rho_633 = v6(Psi3)/v6(Psi6)
    """
    v6_array = vn_array[:, 5]
    nev = len(v6_array)
    v6_Psi6 = sqrt(mean(abs(v6_array)**2.))
    v6_Psi6_err = std(abs(v6_array)**2.)/(2.*v6_Psi6)/sqrt(nev)
    rho_633 = v6_Psi3/v6_Psi6
    rho_633_err = sqrt((v6_Psi3_err/v6_Psi6)**2.
                        + (v6_Psi3*v6_Psi6_err/v6_Psi6**2.)**2.)
    return(rho_633, rho_633_err)


def calculate_v4_L(v4_Psi2, v4_Psi2_err, vn_array):
    """
        v4_L = sqrt(v4(Psi4)^2 - v4(Psi2)^2)
    """
    v4_array = vn_array[:, 3]
    nev = len(v4_array)
    v4_Psi4_sq = mean(abs(v4_array)**2.)
    v4_Psi4_sq_err = std(abs(v4_array)**2.)/sqrt(nev)
    v4_L = sqrt(v4_Psi4_sq - v4_Psi2**2.)
    v4_L_err = (sqrt(v4_Psi4_sq_err**2. + (2.*v4_Psi2*v4_Psi2_err)**2.)
                /(2.*v4_L))
    return(v4_L, v4_L_err)


def calculate_v5_L(v5_Psi23, v5_Psi23_err, vn_array):
    """
        v5_L = sqrt(v5(Psi5)^2 - v5(Psi23)^2)
    """
    v5_array = vn_array[:, 4]
    nev = len(v5_array)
    v5_Psi5_sq = mean(abs(v5_array)**2.)
    v5_Psi5_sq_err = std(abs(v5_array)**2.)/sqrt(nev)
    v5_L = sqrt(v5_Psi5_sq - v5_Psi23**2.)
    v5_L_err = (sqrt(v5_Psi5_sq_err**2. + (2.*v5_Psi23*v5_Psi23_err)**2.)
                /(2.*v5_L))
    return(v5_L, v5_L_err)


def calculate_v6_L(chi_6222, chi_6222_err, chi_633, chi_633_err, vn_array):
    """
        v6_L = sqrt(v6(Psi6)^2 - chi_6222^2 v2^6
                    - chi_633^2 v3^4 - 2 Re(chi_6222*chi_633*v2^3 v3^{2*}))
    """
    v6_array = vn_array[:, 5]
    v2_array = vn_array[:, 1]
    v3_array = vn_array[:, 2]
    nev = len(v6_array)
    v6_Psi6_sq = mean(abs(v6_array)**2.)
    v6_Psi6_sq_err = std(abs(v6_array)**2.)/sqrt(nev)
    v2_6 = mean(abs(v2_array)**6.)
    v2_6_err = std(abs(v2_array)**6.)/sqrt(nev)
    v3_4 = mean(abs(v3_array)**4.)
    v3_4_err = std(abs(v3_array)**4.)/sqrt(nev)
    v23 = real(mean(v2_array**3.*conj(v3_array)**2.))
    v23_err = real(std(v2_array**3.*conj(v3_array)**2.))/sqrt(nev)
    v6_L = (v6_Psi6_sq - chi_6222**2.*v2_6 - chi_633**2.*v3_4
            - 2.*chi_6222*chi_633*v23)
    v6_L_err = sqrt(
            v6_Psi6_sq_err**2.
            + (2.*chi_6222*chi_6222_err*v2_6)**2. + (chi_6222**2.*v2_6_err)**2.
            + (2.*chi_633*chi_633_err*v3_4)**2. + (chi_633**2.*v3_4_err)**2.
            + (2.*chi_6222_err*chi_633*v23)**2.
            + (2.*chi_6222*chi_633_err*v23)**2.
            + (2.*chi_6222*chi_633*v23_err)**2.)
    return(v6_L, v6_L_err)


def calculate_nonlinear_reponse(vn_array):
    """
        this function computes all the nonlinear response coefficients
        proposed in the paper arXiv: 1502.02502 up to v6
    """
    chi_422, chi_422_err = calculate_chi_422(vn_array)
    v4_Psi2, v4_Psi2_err = calculate_v4_Psi2(chi_422, chi_422_err, vn_array)
    rho_422, rho_422_err = calculate_rho_422(v4_Psi2, v4_Psi2_err, vn_array)
    v4_L, v4_L_err = calculate_v4_L(v4_Psi2, v4_Psi2_err, vn_array)

    chi_523, chi_523_err = calculate_chi_523(vn_array)
    v5_Psi23, v5_Psi23_err = calculate_v5_Psi23(chi_523, chi_523_err, vn_array)
    rho_523, rho_523_err = calculate_rho_523(v5_Psi23, v5_Psi23_err, vn_array)
    v5_L, v5_L_err = calculate_v5_L(v5_Psi23, v5_Psi23_err, vn_array)

    chi_6222, chi_6222_err = calculate_chi_6222(vn_array)
    v6_Psi2, v6_Psi2_err = calculate_v6_Psi2(chi_6222, chi_6222_err, vn_array)
    rho_6222, rho_6222_err = calculate_rho_6222(v6_Psi2, v6_Psi2_err, vn_array)

    chi_633, chi_633_err = calculate_chi_633(vn_array)
    v6_Psi3, v6_Psi3_err = calculate_v6_Psi3(chi_633, chi_633_err, vn_array)
    rho_633, rho_633_err = calculate_rho_633(v6_Psi3, v6_Psi3_err, vn_array)
    v6_L, v6_L_err = calculate_v6_L(chi_6222, chi_6222_err,
                                    chi_633, chi_633_err, vn_array)

    results = [v4_L, v4_L_err, v4_Psi2, v4_Psi2_err, rho_422, rho_422_err,
               chi_422, chi_422_err,
               v5_L, v5_L_err, v5_Psi23, v5_Psi23_err, rho_523, rho_523_err,
               chi_523, chi_523_err,
               v6_L, v6_L_err, v6_Psi2, v6_Psi2_err, v6_Psi3, v6_Psi3_err,
               rho_6222, rho_6222_err, rho_633, rho_633_err,
               chi_6222, chi_6222_err, chi_633, chi_633_err]
    return(results)


def calcualte_vn_2(vn_data_array):
    vn_data_array = array(vn_data_array)
    nev = len(vn_data_array[:, 0])
    vn_2 = sqrt(mean(abs(vn_data_array)**2., 0)) + 1e-30
    vn_2_err = std(abs(vn_data_array)**2., 0)/sqrt(nev)/2./vn_2
    return(vn_2, vn_2_err)


def calculate_diff_vn_single_event(pT_ref_low, pT_ref_high, data):
    npT = 50
    pT_inte_array = linspace(pT_ref_low, pT_ref_high, npT)
    dN_event = data[:, 2]
    pT_event = data[:, 0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    dN_ref = sum(dN_interp*pT_inte_array)
    temp_vn_real_array = []
    temp_vn_imag_array = []
    temp_vn_denorm_array = []
    for iorder in range(1, n_order):
        vn_real_event = data[:, 3*iorder]
        vn_imag_event = data[:, 3*iorder+1]
        vn_real_interp = interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = interp(pT_inte_array, pT_event, vn_imag_event)
        vn_real_inte = (
            sum(vn_real_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_imag_inte = (
            sum(vn_imag_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_ref = vn_real_inte + 1j*vn_imag_inte
        vn_pt = vn_real_event + 1j*vn_imag_event
        numerator_real = real(dN_event*vn_pt*dN_ref*conj(vn_ref))
        numerator_imag = imag(dN_event*vn_pt*dN_ref*conj(vn_ref))
        denorm = dN_event*dN_ref
        temp_vn_real_array.append(numerator_real)
        temp_vn_imag_array.append(numerator_imag)
    temp_vn_denorm_array.append(denorm)
    return(temp_vn_real_array, temp_vn_imag_array, temp_vn_denorm_array)


def get_vn_diff_2PC_from_single_event(data):
    dN_event = data[:, 2]
    temp_vn_real_array = []
    temp_vn_imag_array = []
    temp_vn_denorm_array = []
    for iorder in range(1, n_order):
        vn_real_event = data[:, 3*iorder]
        vn_imag_event = data[:, 3*iorder+1]
        vn_pt = vn_real_event + 1j*vn_imag_event
        numerator_real = real(dN_event*vn_pt)
        numerator_imag = imag(dN_event*vn_pt)
        denorm = dN_event
        temp_vn_real_array.append(numerator_real)
        temp_vn_imag_array.append(numerator_imag)
    temp_vn_denorm_array.append(denorm)
    return(temp_vn_real_array, temp_vn_imag_array, temp_vn_denorm_array)


def calculate_vn_diff_SP(vn_diff_real, vn_diff_imag, vn_diff_denorm,
                         vn_2, vn_2_err):
    """
        this funciton calculates the scalar-product vn
    """
    vn_diff_real = array(vn_diff_real)
    vn_diff_imag = array(vn_diff_imag)
    vn_diff_denorm = array(vn_diff_denorm) + 1e-30
    nev = len(vn_diff_denorm[:, 0])
    vn_denorm = vn_2.reshape(len(vn_2), 1)
    vn_denorm_err = vn_2_err.reshape(len(vn_2_err), 1)
    vn_diff_SP = (
        mean(vn_diff_real, 0)/mean(vn_diff_denorm, 0)/vn_denorm)
    vn_diff_SP_err = sqrt(
        ( std(vn_diff_real, 0)/sqrt(nev)/mean(vn_diff_denorm, 0)
          /vn_denorm)**2.
        + (vn_diff_SP*vn_denorm_err/vn_denorm)**2.)
    return(vn_diff_SP, vn_diff_SP_err)


def calculate_vn_diff_2PC(vn_diff_real, vn_diff_imag, vn_diff_denorm):
    """
        this funciton calculates the rms vn[2](pT)
    """
    vn_diff_real = array(vn_diff_real)
    vn_diff_imag = array(vn_diff_imag)
    nev = len(vn_diff_real[:, 0])
    vn_diff_2PC = sqrt(mean((vn_diff_real**2. + vn_diff_imag**2.), 0))
    vn_diff_2PC_err = (std((vn_diff_real**2. + vn_diff_imag**2.), 0)/sqrt(nev)
                       /(2.*vn_diff_2PC + 1e-15))
    return(vn_diff_2PC, vn_diff_2PC_err)


def calculate_vn_distribution(vn_array):
    nbin = 20
    vn_array = array(vn_array)
    vn_dim = len(vn_array[0, :])
    output = []
    for vn_order in range(vn_dim):
        vn_mag_array = abs(vn_array[:, vn_order])
        vn_min = min(vn_mag_array)
        vn_max = max(vn_mag_array) + 1e-10
        bin_boundaries = linspace(vn_min, vn_max, nbin+1)
        bin_width = bin_boundaries[1] - bin_boundaries[0]
        bin_center = zeros([nbin])
        bin_value = zeros([nbin])
        for vn_elem in vn_mag_array:
            vn_idx = int((vn_elem - vn_min)/bin_width)
            bin_value[vn_idx] += 1.
            bin_center[vn_idx] += vn_elem
        bin_center = bin_center/(bin_value + 1e-15)
        bin_value = bin_value/len(vn_array)
        bin_value_err = sqrt(bin_value/len(vn_array))
        bin_value = bin_value/bin_width
        bin_value_err = bin_value_err/bin_width
        for i in range(nbin):
            if abs(bin_center[i]) < 1e-15:
                bin_center[i] = (bin_boundaries[i] + bin_boundaries[i+1])/2.
        output.append(bin_center)
        output.append(bin_value)
        output.append(bin_value_err)
    output = array(output)
    return(output.transpose())


def calcualte_event_plane_correlations(vn_array):
    """
        this function compute the scalar-product event plane correlations
        vn_array is a matrix [event_idx, vn_order]
    """
    vn_array = array(vn_array)
    nev = len(vn_array[:, 0])
    # cos(4(Psi_2 - Psi_4))
    v2_array = vn_array[:, 1]
    v4_array = vn_array[:, 3]
    v2_2 = mean(abs(v2_array)**2.)
    v4_2 = mean(abs(v4_array)**2.)
    v2_2_err = std(abs(v2_array)**2.)/sqrt(nev)
    v4_2_err = std(abs(v4_array)**2.)/sqrt(nev)
    corr_224_num = mean(real(v2_array**2.*conj(v4_array)))
    corr_224_num_err = std(real(v2_array**2.*conj(v4_array)))/sqrt(nev)
    corr_224_denorm = sqrt(v2_2**2.*v4_2)
    corr_224_denorm_err = sqrt((v2_2_err/(v2_2**2.*sqrt(v4_2)))**2.
                               + (v4_2_err/(2.*v2_2*v4_2**1.5))**2.)
    corr_224 = corr_224_num/sqrt(v2_2*v2_2*v4_2)
    corr_224_err = sqrt(
        (corr_224_num_err/corr_224_denorm)**2.
         + (corr_224_denorm_err*corr_224_num)**2.)

    # cos(6(Psi_2 - Psi_3))
    v3_array = vn_array[:, 2]
    v3_2 = mean(abs(v3_array)**2.)
    v3_2_err = std(abs(v3_array)**2.)/sqrt(nev)
    corr_22233_num = mean(real(v2_array**3.*conj(v3_array)**2.))
    corr_22233_num_err = std(real(v2_array**3.*conj(v3_array)**2.))/sqrt(nev)
    corr_22233_denorm = sqrt(v2_2**3.*v3_2**2.)
    corr_22233_denorm_err = sqrt((3.*v2_2_err/(2.*v2_2**2.5*v3_2))**2.
                                 + (v3_2_err/(v2_2**1.5*v3_2**2.))**2.)
    corr_22233 = corr_22233_num/corr_22233_denorm
    corr_22233_err = sqrt(
        (corr_22233_num_err/corr_22233_denorm)**2.
        + (corr_22233_denorm_err*corr_22233_num)**2.)

    # cos(6(Psi_2 - Psi_6))
    v6_array = vn_array[:, 5]
    v6_2 = mean(abs(v6_array)**2.)
    v6_2_err = std(abs(v6_array)**2.)/sqrt(nev)
    corr_2226_num = mean(real(v2_array**3.*conj(v6_array)))
    corr_2226_num_err = std(real(v2_array**3.*conj(v6_array)))/sqrt(nev)
    corr_2226_denorm = sqrt(v2_2**3.*v6_2)
    corr_2226_denorm_err = sqrt(
        (3.*v2_2_err/(2.*v2_2**2.5*sqrt(v6_2)))**2.
        + (v6_2_err/(2.*v2_2**1.5*v6_2**1.5))**2.)
    corr_2226 = corr_2226_num/corr_2226_denorm
    corr_2226_err = sqrt(
        (corr_2226_num_err/corr_2226_denorm)**2.
        + (corr_2226_num*corr_2226_denorm_err)**2.)

    # cos(6(Psi_3 - Psi_6))
    corr_336_num = mean(real(v3_array**2.*conj(v6_array)))
    corr_336_num_err = std(real(v3_array**2.*conj(v6_array)))/sqrt(nev)
    corr_336_denorm = sqrt(v3_2**2.*v6_2)
    corr_336_denorm_err = sqrt(
            (v3_2_err/(v3_2**2.*sqrt(v6_2)))**2.
            + (v6_2_err/(2.*v3_2*v6_2**1.5))**2.)
    corr_336 = corr_336_num/corr_336_denorm
    corr_336_err = sqrt(
            (corr_336_num_err/corr_336_denorm)**2.
            + (corr_336_num*corr_336_denorm_err)**2.)

    # cos(2Psi_2 + 3Psi_3 - 5Psi_5)
    v5_array = vn_array[:, 4]
    v5_2 = mean(abs(v5_array)**2.)
    v5_2_err = std(abs(v5_array)**2.)/sqrt(nev)
    corr_235_num = mean(real(v2_array*v3_array*conj(v5_array)))
    corr_235_num_err = std(real(v2_array*v3_array*conj(v5_array)))/sqrt(nev)
    corr_235_denorm = sqrt(v2_2*v3_2*v5_2)
    corr_235_denorm_err = sqrt((v2_2_err/(2.*v2_2*sqrt(v3_2*v5_2)))**2.
                               + (v3_2_err/(2.*v3_2*sqrt(v2_2*v5_2)))**2.
                               + (v5_2_err/(2.*v5_2*sqrt(v2_2*v3_2)))**2.)
    corr_235 = corr_235_num/corr_235_denorm
    corr_235_err = sqrt(
            (corr_235_num_err/corr_235_denorm)**2.
            + (corr_235_num*corr_235_denorm_err)**2.)

    # cos(2Psi_2 + 4Psi_4 - 6Psi_6)
    corr_246_num = mean(real(v2_array*v4_array*conj(v6_array)))
    corr_246_num_err = std(real(v2_array*v4_array*conj(v6_array)))/sqrt(nev)
    corr_246_denorm = sqrt(v2_2*v4_2*v6_2)
    corr_246_denorm_err = sqrt((v2_2_err/(2.*v2_2*sqrt(v4_2*v6_2)))**2.
                               + (v4_2_err/(2.*v4_2*sqrt(v2_2*v6_2)))**2.
                               + (v6_2_err/(2.*v6_2*sqrt(v2_2*v4_2)))**2.)
    corr_246 = corr_246_num/corr_246_denorm
    corr_246_err = sqrt(
        (corr_246_num_err/corr_246_denorm)**2.
        + (corr_246_num*corr_246_denorm_err)**2.)

    # cos(2Psi_2 - 6Psi_3 + 4Psi_4)
    corr_234_num = mean(real(v2_array*conj(v3_array)**2.*v4_array))
    corr_234_num_err = std(real(v2_array*conj(v3_array)**2.*v4_array))/sqrt(nev)
    corr_234_denorm = sqrt(v2_2*v3_2**2.*v4_2)
    corr_234_denorm_err = sqrt(
            (v2_2_err/(2.*v2_2**1.5*v3_2*sqrt(v4_2)))**2.
            + (v3_2_err/(sqrt(v2_2*v4_2)*v3_2**2.))**2.
            + (v4_2_err/(2.*sqrt(v2_2)*v3_2*v4_2**1.5))**2.)
    corr_234 = corr_234_num/corr_234_denorm
    corr_234_err = sqrt(
        (corr_234_num_err/corr_234_denorm)**2.
        + (corr_234_num*corr_234_denorm_err)**2.)

    results = [corr_224, corr_22233, corr_2226, corr_336,
               corr_235, corr_246, corr_234]
    results_err = [corr_224_err, corr_22233_err, corr_2226_err, corr_336_err,
                   corr_235_err, corr_246_err, corr_234_err]
    return(results, results_err)


def calculate_vn_arrays_for_rn_ratios(data):
    # this function compute the complex pT-integrated Vn vector
    # in different pT ranges for a single event
    # it returns a 2d matrix vn_arrays[pT_idx, n_order_idx]
    pT_boundaries = [0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0]
    npT = 50
    vn_arrays = []
    for ipT in range(len(pT_boundaries)-1):
        pT_low = pT_boundaries[ipT]
        pT_high = pT_boundaries[ipT + 1]
        pT_mid = (pT_low + pT_high)/2.
        vn_array = calcualte_inte_vn(pT_low, pT_high, data)
        vn_array.insert(0, pT_mid)
        vn_arrays.append(vn_array)
    return(vn_arrays)


def calculate_rn_ratios(vn_event_arrays):
    # this function compute rn ratio in different pT bins
    # according to the CMS measurements
    # it reads in a 3d data cube
    #       vn_event_arrays[event_idx, pT_idx, n_order_idx]
    # it returns rn_arrays[iorder, pT_idx, 3]
    vn_event_arrays = array(vn_event_arrays)
    rn_arrays = []
    for iorder in range(2, 5):
        # compute r2, r3, r4
        rn_array = []
        for itrig in range(3, len(vn_event_arrays[0, :, 0])):
            pT_trig = real(vn_event_arrays[0, itrig, 0])
            vn_trig_array = vn_event_arrays[:, itrig, iorder]
            nev = len(vn_trig_array)
            denorm2_array = abs(vn_trig_array)**2.
            denorm2 = sqrt(mean(denorm2_array))
            denorm2_err = std(denorm2_array)/sqrt(nev)/(2.*denorm2)
            for iasso in range(0, itrig+1):
                pT_asso = real(vn_event_arrays[0, iasso, 0])
                vn_asso_array = vn_event_arrays[:, iasso, iorder]
                num_array = real(vn_asso_array*conj(vn_trig_array))
                num = mean(num_array)
                num_err = std(num_array)/sqrt(nev)
                denorm1_array = abs(vn_asso_array)**2.
                denorm1 = sqrt(mean(denorm1_array))
                denorm1_err = std(denorm1_array)/sqrt(nev)/(2.*denorm1)

                rn_temp = num/(denorm1*denorm2)
                rn_temp_err = sqrt(
                    (num_err/(denorm1*denorm2))**2.
                    + (num*denorm1_err/((denorm1**2.)*denorm2))**2.
                    + (num*denorm2_err/(denorm1*(denorm2**2.)))**2.)
                rn_array.append([pT_trig - pT_asso, rn_temp, rn_temp_err])
        rn_arrays.append(rn_array)
    rn_arrays = array(rn_arrays)
    return(rn_arrays)


file_folder_list = glob(path.join(working_folder, '*'))
nev = len(file_folder_list)
for ipart, particle_id in enumerate(particle_list):
    print("processing %s ..." % particle_list[ipart])
    
    # first particle yield dN/dy
    file_name = '%s_integrated_vndata.dat' % particle_id

    dN_dy = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))

        dN_dy.append(temp_data[0, 1])

    dN_dy = array(dN_dy)
    dN_dy_avg = mean(dN_dy)
    dN_dy_avg_err = std(dN_dy)/sqrt(nev)

    # then <pT>, vn, dN/(2pi dy pT dpT), vn{SP}(pT)
    file_name = '%s_vndata.dat' % particle_id
   
    pT_array = []
    dN_array = []
    vn_phenix_array = []
    vn_star_array = []
    vn_alice_array = []
    vn_cms_array = []
    vn_cms_arrays_for_rn = []
    vn_atlas_array = []
    vn_diff_phenix_real = []; vn_diff_phenix_imag = [];
    vn_diff_phenix_denorm = []
    vn_diff_star_real = []; vn_diff_star_imag = []; vn_diff_star_denorm = []
    vn_diff_alice_real = []; vn_diff_alice_imag = []; vn_diff_alice_denorm = []
    vn_diff_2PC_real = []; vn_diff_2PC_imag = []; vn_diff_2PC_denorm = []
    vn_diff_cms_real = []; vn_diff_cms_imag = []; vn_diff_cms_denorm = []
    vn_diff_atlas_real = []; vn_diff_atlas_imag = []; vn_diff_atlas_denorm = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))
        
        dN_event = temp_data[:, 2]  # dN/(2pi dy pT dpT)
        pT_event = temp_data[:, 0]

        # record particle spectra
        pT_array.append(pT_event)
        dN_array.append(dN_event)

        # pT-integrated vn
        # vn with PHENIX pT cut
        temp_vn_array = calcualte_inte_vn(0.2, 2.0, temp_data)
        vn_phenix_array.append(temp_vn_array)

        # vn with STAR pT cut
        temp_vn_array = calcualte_inte_vn(0.15, 2.0, temp_data)
        vn_star_array.append(temp_vn_array)

        # vn with ALICE pT cut
        temp_vn_array = calcualte_inte_vn(0.2, 3.0, temp_data)
        vn_alice_array.append(temp_vn_array)
        
        # vn with CMS pT cut
        temp_vn_array = calcualte_inte_vn(0.3, 3.0, temp_data)
        vn_cms_array.append(temp_vn_array)
        if particle_id == "Charged_eta":
            temp_vn_arrays = (
                    calculate_vn_arrays_for_rn_ratios(temp_data))
            vn_cms_arrays_for_rn.append(temp_vn_arrays)
        
        # vn with ATLAS pT cut
        temp_vn_array = calcualte_inte_vn(0.5, 3.0, temp_data)
        vn_atlas_array.append(temp_vn_array)

        # pT-differential vn using scalar-product method
        # vn{SP}(pT) with PHENIX pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                        calculate_diff_vn_single_event(0.15, 2.0, temp_data))
        vn_diff_phenix_real.append(temp_vn_diff_real);
        vn_diff_phenix_imag.append(temp_vn_diff_imag);
        vn_diff_phenix_denorm.append(temp_dn_diff);

        # vn{SP}(pT) with STAR pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                        calculate_diff_vn_single_event(0.15, 2.0, temp_data))
        vn_diff_star_real.append(temp_vn_diff_real);
        vn_diff_star_imag.append(temp_vn_diff_imag);
        vn_diff_star_denorm.append(temp_dn_diff);
        
        # vn{SP}(pT) with ALICE pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                        calculate_diff_vn_single_event(0.2, 3.0, temp_data))
        vn_diff_alice_real.append(temp_vn_diff_real);
        vn_diff_alice_imag.append(temp_vn_diff_imag);
        vn_diff_alice_denorm.append(temp_dn_diff);
        
        # vn{SP}(pT) with CMS pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                        calculate_diff_vn_single_event(0.3, 3.0, temp_data))
        vn_diff_cms_real.append(temp_vn_diff_real);
        vn_diff_cms_imag.append(temp_vn_diff_imag);
        vn_diff_cms_denorm.append(temp_dn_diff);
        
        # vn{SP}(pT) with ATLAS pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                        calculate_diff_vn_single_event(0.5, 3.0, temp_data))
        vn_diff_atlas_real.append(temp_vn_diff_real);
        vn_diff_atlas_imag.append(temp_vn_diff_imag);
        vn_diff_atlas_denorm.append(temp_dn_diff);
        
        # pT-differential vn using 2PC method
        # vn[2](pT)
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                                get_vn_diff_2PC_from_single_event(temp_data))
        vn_diff_2PC_real.append(temp_vn_diff_real)
        vn_diff_2PC_imag.append(temp_vn_diff_imag)
        vn_diff_2PC_denorm.append(temp_dn_diff)

    # now we perform event average
    dN_array = array(dN_array)
    pT_array = array(pT_array)
    n_pT = len(pT_array[0, :])
    pT_spectra = zeros([n_pT])
    for ipT in range(len(pT_array[0, :])):
        dN_temp = sum(dN_array[:, ipT]*pT_array[:, ipT])
        if(dN_temp > 0):
            pT_spectra[ipT] = (
                    sum(pT_array[:, ipT]**2.*dN_array[:, ipT])/dN_temp)
        else:
            pT_spectra[ipT] = mean(pT_array[:, ipT])
    dN_spectra = mean(pT_array*dN_array, 0)/pT_spectra   # dN/(2pi dy pT dpT)
    dN_spectra_err = std(pT_array*dN_array, 0)/pT_spectra/sqrt(nev)

    # calculate mean pT
    pT_interp = linspace(0.05, 2.95, 30)
    dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra+1e-30)))
    dN_interp_err = interp(pT_interp, pT_spectra, dN_spectra_err)
    mean_pT = sum(pT_interp**2.*dN_interp)/sum(pT_interp*dN_interp)
    mean_pT_upper = (sum(pT_interp**2.*(dN_interp+dN_interp_err))
                     /sum(pT_interp*(dN_interp+dN_interp_err)))
    mean_pT_lower = (sum(pT_interp**2.*(dN_interp-dN_interp_err))
                     /sum(pT_interp*(dN_interp-dN_interp_err)))
    mean_pT_err = max(abs(mean_pT_upper - mean_pT), 
                      abs(mean_pT - mean_pT_lower))

    # calcualte vn{2}
    vn_phenix_2, vn_phenix_2_err = calcualte_vn_2(vn_phenix_array)
    vn_star_2, vn_star_2_err = calcualte_vn_2(vn_star_array)
    vn_alice_2, vn_alice_2_err = calcualte_vn_2(vn_alice_array)
    vn_cms_2, vn_cms_2_err = calcualte_vn_2(vn_cms_array)
    vn_atlas_2, vn_atlas_2_err = calcualte_vn_2(vn_atlas_array)

    if (particle_id == 'Charged_eta'):
        vn_alice_array2 = array(vn_alice_array)
        vn_cms_array2 = array(vn_cms_array)
        vn_atlas_array2 = array(vn_atlas_array)
        # calculate non-linear response coefficents with ALICE pT cut
        nonlinear_response_alice = calculate_nonlinear_reponse(vn_alice_array2)
        # calculate non-linear response coefficents with CMS pT cut
        nonlinear_response_cms = calculate_nonlinear_reponse(vn_cms_array2)
        # calculate non-linear response coefficents with ATLAS pT cut
        nonlinear_response_atlas = calculate_nonlinear_reponse(vn_atlas_array2)

        # calculate vn distribution for charged hadrons
        vn_phenix_dis = calculate_vn_distribution(vn_phenix_array)
        vn_star_dis = calculate_vn_distribution(vn_star_array)
        vn_alice_dis = calculate_vn_distribution(vn_alice_array)
        vn_cms_dis = calculate_vn_distribution(vn_cms_array)
        vn_atlas_dis = calculate_vn_distribution(vn_atlas_array)

        # calculate rn ratios
        rn_cms = calculate_rn_ratios(vn_cms_arrays_for_rn)

        # calculate flow event-plane correlation
        vn_corr_alice, vn_corr_alice_err = (
                calcualte_event_plane_correlations(vn_alice_array))
        vn_corr_atlas, vn_corr_atlas_err = (
                calcualte_event_plane_correlations(vn_atlas_array))

    # calcualte vn{SP}(pT)
    vn_diff_SP_phenix, vn_diff_SP_phenix_err = calculate_vn_diff_SP(
            vn_diff_phenix_real, vn_diff_phenix_imag, vn_diff_phenix_denorm,
            vn_phenix_2, vn_phenix_2_err)

    vn_diff_SP_star, vn_diff_SP_star_err = calculate_vn_diff_SP(
            vn_diff_star_real, vn_diff_star_imag, vn_diff_star_denorm,
            vn_star_2, vn_star_2_err)
    
    vn_diff_SP_alice, vn_diff_SP_alice_err = calculate_vn_diff_SP(
            vn_diff_alice_real, vn_diff_alice_imag, vn_diff_alice_denorm,
            vn_alice_2, vn_alice_2_err)
    
    vn_diff_SP_cms, vn_diff_SP_cms_err = calculate_vn_diff_SP(
            vn_diff_cms_real, vn_diff_cms_imag, vn_diff_cms_denorm,
            vn_cms_2, vn_cms_2_err)
        
    vn_diff_SP_atlas, vn_diff_SP_atlas_err = calculate_vn_diff_SP(
            vn_diff_atlas_real, vn_diff_atlas_imag, vn_diff_atlas_denorm,
            vn_atlas_2, vn_atlas_2_err)
    
    # calcualte vn[2](pT)
    vn_diff_2PC, vn_diff_2PC_err = calculate_vn_diff_2PC(
            vn_diff_2PC_real, vn_diff_2PC_imag, vn_diff_2PC_denorm)
    
    ###########################################################################
    # finally, output all the results
    ###########################################################################
    
    if (particle_id =='Charged_eta'):
        # output non-linear response coefficients chi_n for CMS pt cut
        output_filename = ("non_linear_response_coefficients_CMS.dat")
        f = open(output_filename, 'w')
        f.write("# type  value  stat. err\n")
        for i in range(len(nonlinear_reponse_correlator_name_list)):
            f.write("%s  %.10e  %.10e\n"
                    % (nonlinear_reponse_correlator_name_list[i],
                       nonlinear_response_cms[2*i],
                       nonlinear_response_cms[2*i+1]))
        f.close()
        shutil.move(output_filename, avg_folder)

        # output non-linear response coefficients chi_n for ALICE pt cut
        output_filename = ("non_linear_response_coefficients_ALICE.dat")
        f = open(output_filename, 'w')
        f.write("# type  value  stat. err\n")
        for i in range(len(nonlinear_reponse_correlator_name_list)):
            f.write("%s  %.10e  %.10e\n"
                    % (nonlinear_reponse_correlator_name_list[i],
                       nonlinear_response_alice[2*i],
                       nonlinear_response_alice[2*i+1]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output non-linear response coefficients chi_n for ATLAS pt cut
        output_filename = ("non_linear_response_coefficients_ATLAS.dat")
        f = open(output_filename, 'w')
        f.write("# type  value  stat. err\n")
        for i in range(len(nonlinear_reponse_correlator_name_list)):
            f.write("%s  %.10e  %.10e\n"
                    % (nonlinear_reponse_correlator_name_list[i],
                       nonlinear_response_atlas[2*i],
                       nonlinear_response_atlas[2*i+1]))
        f.close()
        shutil.move(output_filename, avg_folder)

    output_filename = ("%s_integrated_observables.dat"
                       % particle_list[ipart])
    f = open(output_filename, 'w')
    f.write("dN/dy= %.10e +/- %.10e\n" % (dN_dy_avg, dN_dy_avg_err))
    f.write("<pT>= %.10e +/- %.10e\n" % (mean_pT, mean_pT_err))
    for iorder in range(1, n_order):
        f.write("v_%d{2}(phenix)= %.10e +/- %.10e\n" 
                % (iorder, vn_phenix_2[iorder-1], vn_phenix_2_err[iorder-1]))
        f.write("v_%d{2}(STAR)= %.10e +/- %.10e\n" 
                % (iorder, vn_star_2[iorder-1], vn_star_2_err[iorder-1]))
        f.write("v_%d{2}(ALICE)= %.10e +/- %.10e\n" 
                % (iorder, vn_alice_2[iorder-1], vn_alice_2_err[iorder-1]))
        f.write("v_%d{2}(CMS)= %.10e +/- %.10e\n" 
                % (iorder, vn_cms_2[iorder-1], vn_cms_2_err[iorder-1]))
        f.write("v_%d{2}(ATLAS)= %.10e +/- %.10e\n" 
                % (iorder, vn_atlas_2[iorder-1], vn_atlas_2_err[iorder-1]))
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_PHENIX.dat" 
                       % particle_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_phenix[iorder-1, ipT], 
                                        vn_diff_SP_phenix_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)

    output_filename = ("%s_differential_observables_STAR.dat" 
                       % particle_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_star[iorder-1, ipT], 
                                        vn_diff_SP_star_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_ALICE.dat" 
                       % particle_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_alice[iorder-1, ipT], 
                                        vn_diff_SP_alice_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_2PC.dat" 
                       % particle_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            "vn[2]  vn[2]_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_2PC[iorder-1, ipT], 
                                        vn_diff_2PC_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_CMS.dat" 
                       % particle_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_cms[iorder-1, ipT], 
                                        vn_diff_SP_cms_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_ATLAS.dat" 
                       % particle_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  "
            "vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_atlas[iorder-1, ipT], 
                                        vn_diff_SP_atlas_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)

    if (particle_id == 'Charged_eta'):
        output_filename = ("%s_vn_distribution_PHENIX.dat"
                           % particle_list[ipart])
        f = open(output_filename, 'w')
        f.write("#vn  dP(vn)/dvn  dP(vn)/dvn_err\n")
        for ipT in range(len(vn_phenix_dis[:, 0])):
            for iorder in range(1, n_order):
                f.write("%.10e  %.10e  %.10e  "
                        % (vn_phenix_dis[ipT, 3*(iorder-1)], 
                           vn_phenix_dis[ipT, 3*(iorder-1)+1],
                           vn_phenix_dis[ipT, 3*(iorder-1)+2]))
            f.write("\n")
        f.close()
        shutil.move(output_filename, avg_folder)
        
        output_filename = ("%s_vn_distribution_STAR.dat"
                           % particle_list[ipart])
        f = open(output_filename, 'w')
        f.write("#vn  dP(vn)/dvn  dP(vn)/dvn_err\n")
        for ipT in range(len(vn_star_dis[:, 0])):
            for iorder in range(1, n_order):
                f.write("%.10e  %.10e  %.10e  "
                        % (vn_star_dis[ipT, 3*(iorder-1)], 
                           vn_star_dis[ipT, 3*(iorder-1)+1],
                           vn_star_dis[ipT, 3*(iorder-1)+2]))
            f.write("\n")
        f.close()
        shutil.move(output_filename, avg_folder)
        
        output_filename = ("%s_vn_distribution_ALICE.dat"
                           % particle_list[ipart])
        f = open(output_filename, 'w')
        f.write("#vn  dP(vn)/dvn  dP(vn)/dvn_err\n")
        for ipT in range(len(vn_alice_dis[:, 0])):
            for iorder in range(1, n_order):
                f.write("%.10e  %.10e  %.10e  "
                        % (vn_alice_dis[ipT, 3*(iorder-1)], 
                           vn_alice_dis[ipT, 3*(iorder-1)+1],
                           vn_alice_dis[ipT, 3*(iorder-1)+2]))
            f.write("\n")
        f.close()
        shutil.move(output_filename, avg_folder)
        
        output_filename = ("%s_vn_distribution_CMS.dat"
                           % particle_list[ipart])
        f = open(output_filename, 'w')
        f.write("#vn  dP(vn)/dvn  dP(vn)/dvn_err\n")
        for ipT in range(len(vn_cms_dis[:, 0])):
            for iorder in range(1, n_order):
                f.write("%.10e  %.10e  %.10e  "
                        % (vn_cms_dis[ipT, 3*(iorder-1)], 
                           vn_cms_dis[ipT, 3*(iorder-1)+1],
                           vn_cms_dis[ipT, 3*(iorder-1)+2]))
            f.write("\n")
        f.close()
        shutil.move(output_filename, avg_folder)
        
        output_filename = ("%s_vn_distribution_ATLAS.dat"
                           % particle_list[ipart])
        f = open(output_filename, 'w')
        f.write("#vn  dP(vn)/dvn  dP(vn)/dvn_err\n")
        for ipT in range(len(vn_atlas_dis[:, 0])):
            for iorder in range(1, n_order):
                f.write("%.10e  %.10e  %.10e  "
                        % (vn_atlas_dis[ipT, 3*(iorder-1)], 
                           vn_atlas_dis[ipT, 3*(iorder-1)+1],
                           vn_atlas_dis[ipT, 3*(iorder-1)+2]))
            f.write("\n")
        f.close()
        shutil.move(output_filename, avg_folder)

        # output rn ratios
        pT_trig = ['1.0', '1.5', '2.0', '2.5', '3.0']
        ipTtrig = 0
        output_filename = ("%s_rn_ratios_CMS_pTtrig_%s_%s.dat"
                           % (particle_list[ipart],
                              pT_trig[ipTtrig], pT_trig[ipTtrig+1]))
        f = open(output_filename, 'w')
        f.write("#pT_mid  rn  rn_err (n = 2, 3, 4)\n")
        for ipT in range(len(rn_cms[0, :, 0])):
            for iorder in range(len(rn_cms[:, 0, 0])):
                f.write("%.5e  %.5e  %.5e  "
                        % (rn_cms[iorder, ipT, 0],
                           rn_cms[iorder, ipT, 1],
                           rn_cms[iorder, ipT, 2]))
            f.write("\n")
            if rn_cms[0, ipT, 0] == 0.0:
                f.close()
                shutil.move(output_filename, avg_folder)
                ipTtrig += 1
                if ipTtrig < (len(pT_trig) - 1):
                    output_filename = ("%s_rn_ratios_CMS_pTtrig_%s_%s.dat"
                                       % (particle_list[ipart],
                                          pT_trig[ipTtrig],
                                          pT_trig[ipTtrig+1]))
                    f = open(output_filename, 'w')
                    f.write("#pT_mid  rn  rn_err (n = 2, 3, 4)\n")
        
        # output flow event-plane correlation
        output_filename = ("%s_event_plane_correlation_ALICE.dat"
                           % particle_list[ipart])
        f = open(output_filename, 'w')
        f.write("#correlator  value  value_err\n")
        f.write("4(24)  %.5e  %.5e\n"
                % (vn_corr_alice[0], vn_corr_alice_err[0]))
        f.write("6(23)  %.5e  %.5e\n"
                % (vn_corr_alice[1], vn_corr_alice_err[1]))
        f.write("6(26)  %.5e  %.5e\n"
                % (vn_corr_alice[2], vn_corr_alice_err[2]))
        f.write("6(36)  %.5e  %.5e\n"
                % (vn_corr_alice[3], vn_corr_alice_err[3]))
        f.write("(235)  %.5e  %.5e\n"
                % (vn_corr_alice[4], vn_corr_alice_err[4]))
        f.write("(246)  %.5e  %.5e\n"
                % (vn_corr_alice[5], vn_corr_alice_err[5]))
        f.write("(234)  %.5e  %.5e\n"
                % (vn_corr_alice[6], vn_corr_alice_err[6]))
        f.close()
        shutil.move(output_filename, avg_folder)
        output_filename = ("%s_event_plane_correlation_ATLAS.dat"
                           % particle_list[ipart])
        f = open(output_filename, 'w')
        f.write("#correlator  value  value_err\n")
        f.write("4(24)  %.5e  %.5e\n"
                % (vn_corr_atlas[0], vn_corr_atlas_err[0]))
        f.write("6(23)  %.5e  %.5e\n"
                % (vn_corr_atlas[1], vn_corr_atlas_err[1]))
        f.write("6(26)  %.5e  %.5e\n"
                % (vn_corr_atlas[2], vn_corr_atlas_err[2]))
        f.write("6(36)  %.5e  %.5e\n"
                % (vn_corr_atlas[3], vn_corr_atlas_err[3]))
        f.write("(235)  %.5e  %.5e\n"
                % (vn_corr_atlas[4], vn_corr_atlas_err[4]))
        f.write("(246)  %.5e  %.5e\n"
                % (vn_corr_atlas[5], vn_corr_atlas_err[5]))
        f.write("(234)  %.5e  %.5e\n"
                % (vn_corr_atlas[6], vn_corr_atlas_err[6]))
        f.close()
        shutil.move(output_filename, avg_folder)

print "Analysis is done."

