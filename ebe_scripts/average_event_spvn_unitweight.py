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
        var = input("do you want to delete it? [y/N]")
        if 'y' in var:
            shutil.rmtree(avg_folder)
        else:
            print("please choose another folder path~")
            exit(0)
    mkdir(avg_folder)
except IndexError:
    print("Usage: average_event_spvn.py working_folder results_folder")
    exit(1)

particle_list = ['9999', '211', '321', '2212', '-211', '-321', '-2212', 
                 '3122', '-3122', '3312', '-3312', '3334', '-3334',
                 '333']
particle_name_list = ['charged_hadron', 'pion_p', 'kaon_p', 'proton',
                      'pion_m', 'kaon_m', 'anti_proton',
                      'Lambda', 'anti_Lambda', 'Xi_m', 'anti_Xi_p',
                      'Omega', 'anti_Omega', 'phi']
nonlinear_reponse_correlator_name_list = [
                'v4_L', 'v4(Psi2)', 'rho_422', 'chi_422',
                'v5_L', 'v5(Psi23)', 'rho_523', 'chi_523',
                'v6_L', 'v6(Psi2)', 'v6(Psi3)',
                'rho_6222', 'rho_633', 'chi_6222', 'chi_633']
symmetric_cumulant_name_list = ['SC_32', 'SC_42']

n_order = 7

def calcualte_inte_vn(pT_low, pT_high, data):
    """
        this function calculates the pT-integrated vn in a 
        given pT range (pT_low, pT_high) for every event in the data
    """
    npT = 50
    pT_inte_array = linspace(pT_low, pT_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, 2]
    pT_event = data[:, 0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    N_event = data[:, -1]
    N_interp = exp(interp(pT_inte_array, pT_event, log(N_event+1e-30)))
    N = sum(N_interp)*dpT/0.1
    temp_vn_array = [N,]
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


def calculate_chi_422(vn_array):
    """
        chi_422 = Re(v4*conj(v2)**2.)/|v2|^4
    """
    dN = real(vn_array[:, 0])
    Q2 = dN*vn_array[:, 2]
    Q4 = dN*vn_array[:, 4]
    nev = len(dN)
    
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))
    
    N3_weight = dN*(dN - 1.)*(dN - 2.)
    chi_422_num = Q4*conj(Q2)*conj(Q2) - 2.*Q2*conj(Q2) - Q4*conj(Q4) + 2.*dN

    chi_422_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        chi_422_JK[iev] = (
            mean(real(chi_422_num[array_idx])/N3_weight[array_idx])
            /(mean(real(Q2_4[array_idx])/N4_weight[array_idx])))

    chi_422_mean = mean(chi_422_JK)
    chi_422_err = sqrt((nev - 1.)/nev*sum((chi_422_JK - chi_422_mean)**2.))
    return(chi_422_mean, chi_422_err)


def calculate_chi_523(vn_array):
    dN = real(vn_array[:, 0])
    Q1 = dN*vn_array[:, 1]
    Q2 = dN*vn_array[:, 2]
    Q3 = dN*vn_array[:, 3]
    Q5 = dN*vn_array[:, 5]
    nev = len(dN)

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q_32 = ((abs(Q2)**2.)*(abs(Q3)**2.) - 2.*real(Q5*conj(Q2)*conj(Q3))
        - 2.*real(Q3*conj(Q1)*conj(Q2)) + abs(Q5)**2. + abs(Q1)**2.
        - (dN - 4.)*(abs(Q2)**2. + abs(Q3)**2.) + dN*(dN - 6.)
    )
    
    N3_weight = dN*(dN - 1.)*(dN - 2.)
    chi_523_num = (Q5*conj(Q2)*conj(Q3) - Q3*conj(Q3) - Q2*conj(Q2)
                   - Q5*conj(Q5) + 2.*dN)

    chi_523_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        chi_523_JK[iev] = (
            mean(real(chi_523_num[array_idx])/N3_weight[array_idx])
            /(mean(real(Q_32[array_idx])/N4_weight[array_idx])))

    chi_523_mean = mean(chi_523_JK)
    chi_523_err = sqrt((nev - 1.)/nev*sum((chi_523_JK - chi_523_mean)**2.))
    return(chi_523_mean, chi_523_err)


def calculate_chi_6222(vn_array):
    dN = real(vn_array[:, 0])
    Q2 = dN*vn_array[:, 2]
    Q4 = dN*vn_array[:, 4]
    Q6 = dN*vn_array[:, 6]
    nev = len(dN)

    N6_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)*(dN - 4.)*(dN - 5.)
    Q2_6 = (abs(Q2)**6. + 9*(abs(Q4)**2.)*(abs(Q2)**2.)
            - 6.*real(Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2))
            + 4.*real(Q6*conj(Q2)*conj(Q2)*conj(Q2))
            - 12.*real(Q6*conj(Q4)*conj(Q2))
            + 18.*(dN - 4.)*real(Q4*conj(Q2)*conj(Q2))
            + 4.*(abs(Q6)**2.)
            - 9.*(dN - 4.)*((abs(Q2)**4.) + (abs(Q4)**2.))
            + 18.*(dN - 5.)*(dN - 2.)*(abs(Q2)**2.)
            - 6.*dN*(dN - 4.)*(dN - 5.))

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    chi_6222_num = (Q6*conj(Q2)*conj(Q2)*conj(Q2) - 3.*Q6*conj(Q4)*conj(Q2)
                    - 3.*Q4*conj(Q2)*conj(Q2) + 2.*Q6*conj(Q6) + 6.*Q2*conj(Q2)
                    + 3.*Q4*conj(Q4) - 6.*dN)

    chi_6222_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        chi_6222_JK[iev] = (
            mean(real(chi_6222_num[array_idx])/N4_weight[array_idx])
            /(mean(real(Q2_6[array_idx])/N6_weight[array_idx])))

    chi_6222_mean = mean(chi_6222_JK)
    chi_6222_err = sqrt((nev - 1.)/nev*sum((chi_6222_JK - chi_6222_mean)**2.))
    return(chi_6222_mean, chi_6222_err)


def calculate_chi_633(vn_array):
    dN = real(vn_array[:, 0])
    Q3 = dN*vn_array[:, 3]
    Q6 = dN*vn_array[:, 6]
    nev = len(dN)

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q3_4 = ((abs(Q3)**4.) - 2.*real(Q6*conj(Q3)*conj(Q3))
             - 4.*(dN - 2.)*(abs(Q3)**2.) + abs(Q6)**2.
             + 2*dN*(dN - 3.))

    N3_weight = dN*(dN - 1.)*(dN - 2.)
    chi_633_num = Q6*conj(Q3)*conj(Q3) - 2.*Q3*conj(Q3) - Q6*conj(Q6) + 2.*dN
    
    chi_633_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        chi_633_JK[iev] = (
            mean(real(chi_633_num[array_idx])/N3_weight[array_idx])
            /(mean(real(Q3_4[array_idx])/N4_weight[array_idx])))

    chi_633_mean = mean(chi_633_JK)
    chi_633_err = sqrt((nev - 1.)/nev*sum((chi_633_JK - chi_633_mean)**2.))
    return(chi_633_mean, chi_633_err)


def calculate_v4_Psi2(chi_422, chi_422_err, vn_array):
    """
        v4(Psi2) = chi_422*sqrt(<abs(V2)**4>)
    """
    dN = real(vn_array[:, 0])
    Q2 = dN*vn_array[:, 2]
    Q4 = dN*vn_array[:, 4]
    nev = len(dN)

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))

    v2_factor = sqrt(mean(Q2_4/N4_weight))
    v2_factor_err = std(Q2_4/N4_weight)/(2.*v2_factor)/sqrt(nev)

    v4_Psi2 = chi_422*v2_factor
    v4_Psi2_err = sqrt((chi_422_err*v2_factor)**2.
                       + (chi_422*v2_factor_err)**2.)
    return(v4_Psi2, v4_Psi2_err)


def calculate_v5_Psi23(chi_523, chi_523_err, vn_array):
    """
        v5(Psi23) = chi_523*sqrt(<abs(V2)**2*abs(V3)**2>)
    """
    dN = real(vn_array[:, 0])
    Q1 = dN*vn_array[:, 1]
    Q2 = dN*vn_array[:, 2]
    Q3 = dN*vn_array[:, 3]
    Q5 = dN*vn_array[:, 5]
    nev = len(dN)
    
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q_32 = ((abs(Q2)**2.)*(abs(Q3)**2.) - 2.*real(Q5*conj(Q2)*conj(Q3))
        - 2.*real(Q3*conj(Q1)*conj(Q2)) + abs(Q5)**2. + abs(Q1)**2.
        - (dN - 4.)*(abs(Q2)**2. + abs(Q3)**2.) + dN*(dN - 6.)
    )

    v23_factor = sqrt(mean(Q_32/N4_weight))
    v23_factor_err = std(Q_32/N4_weight)/(2.*v23_factor)/sqrt(nev)

    v5_Psi23 = chi_523*v23_factor
    v5_Psi23_err = sqrt((chi_523_err*v23_factor)**2.
                        + (chi_523*v23_factor_err)**2.)
    return(v5_Psi23, v5_Psi23_err)


def calculate_v6_Psi2(chi_6222, chi_6222_err, vn_array):
    """
        v6(Psi2) = chi_6222*sqrt(<abs(V2)**6>)
    """
    dN = real(vn_array[:, 0])
    Q2 = dN*vn_array[:, 2]
    Q4 = dN*vn_array[:, 4]
    Q6 = dN*vn_array[:, 6]
    nev = len(dN)

    N6_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)*(dN - 4.)*(dN - 5.)
    Q2_6 = (abs(Q2)**6. + 9*(abs(Q4)**2.)*(abs(Q2)**2.)
            - 6.*real(Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2))
            + 4.*real(Q6*conj(Q2)*conj(Q2)*conj(Q2))
            - 12.*real(Q6*conj(Q4)*conj(Q2))
            + 18.*(dN - 4.)*real(Q4*conj(Q2)*conj(Q2))
            + 4.*(abs(Q6)**2.)
            - 9.*(dN - 4.)*((abs(Q2)**4.) + (abs(Q4)**2.))
            + 18.*(dN - 5.)*(dN - 2.)*(abs(Q2)**2.)
            - 6.*dN*(dN - 4.)*(dN - 5.))

    v2_factor = sqrt(mean(Q2_6/N6_weight))
    v2_factor_err = std(Q2_6/N6_weight)/(2.*v2_factor)/sqrt(nev)

    v6_Psi2 = chi_6222*v2_factor
    v6_Psi2_err = sqrt((chi_6222_err*v2_factor)**2.
                        + (chi_6222*v2_factor_err)**2.)
    return(v6_Psi2, v6_Psi2_err)


def calculate_v6_Psi3(chi_633, chi_633_err, vn_array):
    """
        v6(Psi3) = chi_633*sqrt(<abs(V3)**4>)
    """
    dN = real(vn_array[:, 0])
    Q3 = dN*vn_array[:, 3]
    Q6 = dN*vn_array[:, 6]
    nev = len(dN)

    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q3_4 = ((abs(Q3)**4.) - 2.*real(Q6*conj(Q3)*conj(Q3))
             - 4.*(dN - 2.)*(abs(Q3)**2.) + abs(Q6)**2.
             + 2*dN*(dN - 3.))

    v3_factor = sqrt(mean(Q3_4/N4_weight))
    v3_factor_err = std(Q3_4/N4_weight)/(2.*v3_factor)/sqrt(nev)

    v6_Psi3 = chi_633*v3_factor
    v6_Psi3_err = sqrt((chi_633_err*v3_factor)**2.
                        + (chi_633*v3_factor_err)**2.)
    return(v6_Psi3, v6_Psi3_err)


def calculate_rho_422(v4_Psi2, v4_Psi2_err, vn_array):
    """
        rho_422 = v4(Psi2)/v4(Psi4)
    """
    dN = real(vn_array[:, 0])
    Q4 = dN*vn_array[:, 4]
    nev = len(dN)

    N2_weight = dN*(dN - 1.)
    Q4_2 = abs(Q4)**2. - dN

    v4_Psi4 = sqrt(mean(Q4_2/N2_weight))
    v4_Psi4_err = std(Q4_2/N2_weight)/(2.*v4_Psi4)/sqrt(nev)

    rho_422 = v4_Psi2/v4_Psi4
    rho_422_err = sqrt((v4_Psi2_err/v4_Psi4)**2.
                       + (v4_Psi2*v4_Psi4_err/v4_Psi4**2.)**2.)
    return(rho_422, rho_422_err)


def calculate_rho_523(v5_Psi23, v5_Psi23_err, vn_array):
    """
        rho_523 = v5(Psi23)/v5(Psi5)
    """
    dN = real(vn_array[:, 0])
    Q5 = dN*vn_array[:, 5]
    nev = len(dN)

    N2_weight = dN*(dN - 1.)
    Q5_2 = abs(Q5)**2. - dN

    v5_Psi5 = sqrt(mean(Q5_2/N2_weight))
    v5_Psi5_err = std(Q5_2/N2_weight)/(2.*v5_Psi5)/sqrt(nev)

    rho_523 = v5_Psi23/v5_Psi5
    rho_523_err = sqrt((v5_Psi23_err/v5_Psi5)**2.
                       + (v5_Psi23*v5_Psi5_err/v5_Psi5**2.)**2.)
    return(rho_523, rho_523_err)


def calculate_rho_6222(v6_Psi2, v6_Psi2_err, vn_array):
    """
        rho_6222 = v6(Psi2)/v6(Psi6)
    """
    dN = real(vn_array[:, 0])
    Q6 = dN*vn_array[:, 6]
    nev = len(dN)

    N2_weight = dN*(dN - 1.)
    Q6_2 = abs(Q6)**2. - dN

    v6_Psi6 = sqrt(mean(Q6_2/N2_weight))
    v6_Psi6_err = std(Q6_2/N2_weight)/(2.*v6_Psi6)/sqrt(nev)

    rho_6222 = v6_Psi2/v6_Psi6
    rho_6222_err = sqrt((v6_Psi2_err/v6_Psi6)**2.
                        + (v6_Psi2*v6_Psi6_err/v6_Psi6**2.)**2.)
    return(rho_6222, rho_6222_err)


def calculate_rho_633(v6_Psi3, v6_Psi3_err, vn_array):
    """
        rho_633 = v6(Psi3)/v6(Psi6)
    """
    dN = real(vn_array[:, 0])
    Q6 = dN*vn_array[:, 6]
    nev = len(dN)

    N2_weight = dN*(dN - 1.)
    Q6_2 = abs(Q6)**2. - dN

    v6_Psi6 = sqrt(mean(Q6_2/N2_weight))
    v6_Psi6_err = std(Q6_2/N2_weight)/(2.*v6_Psi6)/sqrt(nev)

    rho_633 = v6_Psi3/v6_Psi6
    rho_633_err = sqrt((v6_Psi3_err/v6_Psi6)**2.
                        + (v6_Psi3*v6_Psi6_err/v6_Psi6**2.)**2.)
    return(rho_633, rho_633_err)


def calculate_v4_L(v4_Psi2, v4_Psi2_err, vn_array):
    """
        v4_L = sqrt(v4(Psi4)^2 - v4(Psi2)^2)
    """
    dN = real(vn_array[:, 0])
    Q4 = dN*vn_array[:, 4]
    nev = len(dN)
    
    N2_weight = dN*(dN - 1.)
    Q4_2 = abs(Q4)**2. - dN

    v4_Psi4_sq = mean(Q4_2/N2_weight)
    v4_Psi4_sq_err = std(Q4_2/N2_weight)/sqrt(nev)

    v4_L = sqrt(v4_Psi4_sq - v4_Psi2**2.)
    v4_L_err = (sqrt(v4_Psi4_sq_err**2. + (2.*v4_Psi2*v4_Psi2_err)**2.)
                /(2.*v4_L))
    return(v4_L, v4_L_err)


def calculate_v5_L(v5_Psi23, v5_Psi23_err, vn_array):
    """
        v5_L = sqrt(v5(Psi5)^2 - v5(Psi23)^2)
    """
    dN = real(vn_array[:, 0])
    Q5 = dN*vn_array[:, 5]
    nev = len(dN)
    
    N2_weight = dN*(dN - 1.)
    Q5_2 = abs(Q5)**2. - dN

    v5_Psi5_sq = mean(Q5_2/N2_weight)
    v5_Psi5_sq_err = std(Q5_2/N2_weight)/sqrt(nev)

    v5_L = sqrt(v5_Psi5_sq - v5_Psi23**2.)
    v5_L_err = (sqrt(v5_Psi5_sq_err**2. + (2.*v5_Psi23*v5_Psi23_err)**2.)
                /(2.*v5_L))
    return(v5_L, v5_L_err)


def calculate_v6_L(chi_6222, chi_6222_err, chi_633, chi_633_err, vn_array):
    """
        v6_L = sqrt(v6(Psi6)^2 - chi_6222^2 v2^6
                    - chi_633^2 v3^4 - 2 Re(chi_6222*chi_633*v2^3 v3^{2*}))
    """
    dN = real(vn_array[:, 0])
    v2_array = vn_array[:, 2]
    v3_array = vn_array[:, 3]
    v6_array = vn_array[:, 6]
    nev = len(dN)
    
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
    """
        this function computes vn{2} and its stat. err.
        self correlation is substracted
    """
    vn_data_array = array(vn_data_array)
    nev = len(vn_data_array[:, 0])
    dN = real(vn_data_array[:, 0])
    dN = dN.reshape(len(dN), 1)
    Qn_array = dN*vn_data_array[:, 1:]
    corr = 1./(dN*(dN - 1.))*(Qn_array*conj(Qn_array) - dN)
    vn_2 = sqrt(real(mean(corr, 0))) + 1e-30
    vn_2_err = std(real(corr), 0)/sqrt(nev)/2./vn_2
    return(nan_to_num(vn_2), nan_to_num(vn_2_err))


def calculate_diff_vn_single_event(pT_ref_low, pT_ref_high, data, data_ref):
    npT = 50
    pT_inte_array = linspace(pT_ref_low, pT_ref_high, npT)
    dpT = pT_inte_array[1] - pT_inte_array[0]
    dN_event = data[:, -1]
    dN_ref_event = data_ref[:, -1]
    pT_ref_event = data_ref[:, 0]
    dN_ref_interp = exp(interp(pT_inte_array, pT_ref_event,
                               log(dN_ref_event + 1e-30)))
    dN_ref = sum(dN_ref_interp)*dpT
    temp_vn_real_array = []
    temp_vn_imag_array = []
    temp_vn_denorm_array1 = []
    temp_vn_denorm_array2 = []
    for iorder in range(1, n_order):
        vn_real_event = data[:, 4*iorder]
        vn_imag_event = data[:, 4*iorder+2]
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
        vn_ref = vn_ref_real_inte + 1j*vn_ref_imag_inte
        vn_pt = vn_real_event + 1j*vn_imag_event
        numerator_real = real(dN_event*vn_pt*dN_ref*conj(vn_ref))
        numerator_imag = imag(dN_event*vn_pt*dN_ref*conj(vn_ref))
        denorm1 = dN_event
        denorm2 = real(dN_ref*vn_ref*dN_ref*conj(vn_ref))
        temp_vn_real_array.append(numerator_real)
        temp_vn_imag_array.append(numerator_imag)
        temp_vn_denorm_array2.append(denorm2)
    temp_vn_denorm_array1.append(denorm1)
    return(temp_vn_real_array, temp_vn_imag_array,
           temp_vn_denorm_array1, temp_vn_denorm_array2)


def get_vn_diff_2PC_from_single_event(data):
    dN_event = data[:, -1]
    temp_vn_real_array = []
    temp_vn_imag_array = []
    temp_vn_denorm_array = []
    for iorder in range(1, n_order):
        vn_real_event = data[:, 4*iorder]
        vn_imag_event = data[:, 4*iorder+2]
        vn_pt = vn_real_event + 1j*vn_imag_event
        numerator_real = real(dN_event*vn_pt)
        numerator_imag = imag(dN_event*vn_pt)
        denorm = dN_event
        temp_vn_real_array.append(numerator_real)
        temp_vn_imag_array.append(numerator_imag)
    temp_vn_denorm_array.append(denorm)
    return(temp_vn_real_array, temp_vn_imag_array, temp_vn_denorm_array)


def calculate_vn_diff_SP(vn_diff_real, vn_diff_imag,
                         vn_diff_denorm1, vn_diff_denorm2):
    """
        this funciton calculates the scalar-product vn
        assumption: no overlap between particles of interest
                    and reference flow Qn vectors
    """
    vn_diff_real = array(vn_diff_real)
    vn_diff_imag = array(vn_diff_imag)
    vn_diff_denorm1 = array(vn_diff_denorm1) + 1e-30
    vn_diff_denorm2 = array(vn_diff_denorm2)
    nev = len(vn_diff_denorm1[:, 0])

    #vn_denorm = vn_2.reshape(len(vn_2), 1)
    #vn_denorm_err = vn_2_err.reshape(len(vn_2_err), 1)
    dn_diff_denorm = mean(vn_diff_denorm1, 0)
    dn_diff_denorm_err = std(vn_diff_denorm1, 0)/sqrt(nev)
    vn_denorm = sqrt(mean(vn_diff_denorm2, 0))
    vn_denorm_err = (std(vn_diff_denorm2, 0)/sqrt(nev)/2./vn_denorm)
    vn_denorm = vn_denorm.reshape(len(vn_denorm), 1)
    vn_denorm_err = vn_denorm_err.reshape(len(vn_denorm), 1)

    vn_diff_SP = mean(vn_diff_real, 0)/dn_diff_denorm/vn_denorm
    vn_diff_SP_err = sqrt(
        (std(vn_diff_real, 0)/sqrt(nev)/dn_diff_denorm/vn_denorm)**2.
        + (vn_diff_SP*dn_diff_denorm_err/dn_diff_denorm)**2.
        + (vn_diff_SP*vn_denorm_err/vn_denorm)**2.)
    return(vn_diff_SP, vn_diff_SP_err)


def calculate_vn_diff_2PC(vn_diff_real, vn_diff_imag, vn_diff_denorm):
    """
        this funciton calculates the rms vn[2](pT)
    """
    vn_diff_real = array(vn_diff_real)
    vn_diff_imag = array(vn_diff_imag)
    vn_diff_denorm = array(vn_diff_denorm)
    nev = len(vn_diff_denorm[:, 0])
    vn_diff_2PC = sqrt(
        mean((vn_diff_real**2. + vn_diff_imag**2. - vn_diff_denorm)
             /(vn_diff_denorm**2. - vn_diff_denorm + 1e-15), 0))
    vn_diff_2PC_err = (
        std((vn_diff_real**2. + vn_diff_imag**2. - vn_diff_denorm)
            /(vn_diff_denorm**2. - vn_diff_denorm + 1e-15), 0)
        /sqrt(nev)/(2.*vn_diff_2PC + 1e-15))
    return(nan_to_num(vn_diff_2PC), nan_to_num(vn_diff_2PC_err))


def calculate_vn_distribution(vn_array):
    nbin = 20
    vn_array = array(vn_array)
    vn_dim = len(vn_array[0, :])
    output = []
    for vn_order in range(vn_dim):
        vn_mag_array = abs(vn_array[:, vn_order])
        vn_min = min(vn_mag_array)
        vn_max = max(vn_mag_array)*1.0001
        bin_boundaries = linspace(vn_min, vn_max, nbin+1)
        bin_width = bin_boundaries[1] - bin_boundaries[0]
        bin_center = zeros([nbin])
        bin_value = zeros([nbin])
        for vn_elem in vn_mag_array:
            vn_idx = int(floor((vn_elem - vn_min)/bin_width))
            if (vn_idx == 20):
                print(vn_elem, vn_min, bin_width)
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
    v2_array = vn_array[:, 2]
    v3_array = vn_array[:, 3]
    v4_array = vn_array[:, 4]
    v5_array = vn_array[:, 5]
    v6_array = vn_array[:, 6]

    corr_224_JK = zeros(nev)
    corr_22233_JK = zeros(nev)
    corr_2226_JK = zeros(nev)
    corr_336_JK = zeros(nev)
    corr_235_JK = zeros(nev)
    corr_246_JK = zeros(nev)
    corr_234_JK = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        v2_2 = mean(abs(v2_array[array_idx])**2.)
        v3_2 = mean(abs(v3_array[array_idx])**2.)
        v4_2 = mean(abs(v4_array[array_idx])**2.)
        v5_2 = mean(abs(v5_array[array_idx])**2.)
        v6_2 = mean(abs(v6_array[array_idx])**2.)

        # cos(4(Psi_2 - Psi_4))
        corr_224_num = mean(real((v2_array[array_idx]**2.)
                                 *conj(v4_array[array_idx])))
        corr_224_JK[iev] = corr_224_num/sqrt(v2_2*v2_2*v4_2)

        # cos(6(Psi_2 - Psi_3))
        corr_22233_num = mean(real((v2_array[array_idx]**3.)
                                   *conj(v3_array[array_idx])**2.))
        corr_22233_JK[iev] = corr_22233_num/sqrt(v2_2**3.*v3_2**2.)

        # cos(6(Psi_2 - Psi_6))
        corr_2226_num = mean(real(v2_array[array_idx]**3.
                                  *conj(v6_array[array_idx])))
        corr_2226_JK[iev] = corr_2226_num/sqrt((v2_2**3.)*v6_2)

        # cos(6(Psi_3 - Psi_6))
        corr_336_num = mean(real((v3_array[array_idx]**2.)
                                 *conj(v6_array[array_idx])))
        corr_336_JK[iev] = corr_336_num/sqrt((v3_2**2.)*v6_2)

        # cos(2Psi_2 + 3Psi_3 - 5Psi_5)
        corr_235_num = mean(real(v2_array[array_idx]*v3_array[array_idx]
                                 *conj(v5_array[array_idx])))
        corr_235_JK[iev] = corr_235_num/sqrt(v2_2*v3_2*v5_2)

        # cos(2Psi_2 + 4Psi_4 - 6Psi_6)
        corr_246_num = mean(real(v2_array[array_idx]*v4_array[array_idx]
                                 *conj(v6_array[array_idx])))
        corr_246_JK[iev] = corr_246_num/sqrt(v2_2*v4_2*v6_2)

        # cos(2Psi_2 - 6Psi_3 + 4Psi_4)
        corr_234_num = mean(real(v2_array[array_idx]
                                 *(conj(v3_array[array_idx])**2.)
                                 *v4_array[array_idx]))
        corr_234_JK[iev] = corr_234_num/sqrt(v2_2*(v3_2**2.)*v4_2)

    corr_224 = mean(corr_224_JK)
    corr_224_err = sqrt((nev - 1.)/nev*sum((corr_224_JK - corr_224)**2.))
    corr_22233 = mean(corr_22233_JK)
    corr_22233_err = sqrt((nev - 1.)/nev*sum((corr_22233_JK - corr_22233)**2.))
    corr_2226 = mean(corr_2226_JK)
    corr_2226_err = sqrt((nev - 1.)/nev*sum((corr_2226_JK - corr_2226)**2.))
    corr_336 = mean(corr_336_JK)
    corr_336_err = sqrt((nev - 1.)/nev*sum((corr_336_JK - corr_336)**2.))
    corr_235 = mean(corr_235_JK)
    corr_235_err = sqrt((nev - 1.)/nev*sum((corr_235_JK - corr_235)**2.))
    corr_246 = mean(corr_246_JK)
    corr_246_err = sqrt((nev - 1.)/nev*sum((corr_246_JK - corr_246)**2.))
    corr_234 = mean(corr_234_JK)
    corr_234_err = sqrt((nev - 1.)/nev*sum((corr_234_JK - corr_234)**2.))

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
    for iorder in range(3, 6):
        # compute r2, r3, r4
        rn_array = []
        for itrig in range(3, len(vn_event_arrays[0, :, 0])):
            pT_trig = real(vn_event_arrays[0, itrig, 0])
            dN_trig = real(vn_event_arrays[:, itrig, 1])
            Qn_trig_array = dN_trig*vn_event_arrays[:, itrig, iorder]
            nev = len(Qn_trig_array)

            denorm2_dN = dN_trig*(dN_trig - 1.)
            denorm2_array = abs(Qn_trig_array)**2. - dN_trig

            for iasso in range(0, itrig+1):
                pT_asso = real(vn_event_arrays[0, iasso, 0])
                dN_asso = real(vn_event_arrays[:, iasso, 1])
                Qn_asso_array = dN_asso*vn_event_arrays[:, iasso, iorder]

                num_dN = dN_trig*dN_asso
                num_array = real(Qn_asso_array*conj(Qn_trig_array))
                if iasso == itrig:
                    num_dN -= dN_asso
                    num_array = (real(Qn_asso_array*conj(Qn_trig_array))
                                 - dN_asso)

                denorm1_dN = dN_asso*(dN_asso - 1.)
                denorm1_array = abs(Qn_asso_array)**2. - dN_asso

                rn_jackknife = zeros(nev)
                for iev in range(nev):
                    array_idx = [True]*nev
                    array_idx[iev] = False
                    array_idx = array(array_idx)

                    num = mean(num_array[array_idx]/num_dN[array_idx])
                    denorm1 = mean(denorm1_array[array_idx]
                                   /denorm1_dN[array_idx])
                    denorm2 = mean(denorm2_array[array_idx]
                                   /denorm2_dN[array_idx])

                    if denorm1 > 0. and denorm2 > 0.:
                        rn_jackknife[iev] = num/sqrt(denorm1*denorm2)

                rn_mean = mean(rn_jackknife)
                rn_err = sqrt((nev - 1.)/nev*sum((rn_jackknife - rn_mean)**2.))
                rn_array.append([pT_trig - pT_asso, rn_mean, rn_err])
        rn_arrays.append(rn_array)
    rn_arrays = array(rn_arrays)
    return(rn_arrays)


def calculate_symmetric_cumulant(vn_data_array):
    """
        this funciton computes the symmetric cumulant
            SC(m,n) = <v_m*conj(v_m)*v_n*conj(v_n)>
                      - <v_m*conj(v_m)>*<v_n*conj(v_n)>
        we use Jackknife resampling method to estimate the statistical error
    """
    vn_data_array = array(vn_data_array)
    nev = len(vn_data_array[:, 0])
    dN = real(vn_data_array[:, 0])
    Q1 = dN*vn_data_array[:, 1]
    Q2 = dN*vn_data_array[:, 2]
    Q3 = dN*vn_data_array[:, 3]
    Q4 = dN*vn_data_array[:, 4]
    Q5 = dN*vn_data_array[:, 5]
    Q6 = dN*vn_data_array[:, 6]

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
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        # SC(3,2)
        SC32_array[iev] = (mean(Q_32[array_idx]/N4_weight[array_idx])
                           - mean(Q3_2[array_idx]/N2_weight[array_idx])
                             *mean(Q2_2[array_idx]/N2_weight[array_idx]))

        # SC(4,2)
        SC42_array[iev] = (mean(Q_42[array_idx]/N4_weight[array_idx])
                           - mean(Q4_2[array_idx]/N2_weight[array_idx])
                             *mean(Q2_2[array_idx]/N2_weight[array_idx]))

    SC32_mean = mean(SC32_array)
    SC32_err = sqrt((nev - 1.)/nev*sum((SC32_array - SC32_mean)**2.))

    SC42_mean = mean(SC42_array)
    SC42_err = sqrt((nev - 1.)/nev*sum((SC42_array - SC42_mean)**2.))

    results = [SC32_mean, SC32_err, SC42_mean, SC42_err]
    return(results)


def calculate_vn4(vn_data_array):
    """
        this funciton computes the 4 particle cumulant vn{4}
            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)
    """
    vn_data_array = array(vn_data_array)
    nev = len(vn_data_array[:, 0])
    dN = real(vn_data_array[:, 0])
    Q1 = dN*vn_data_array[:, 1]
    Q2 = dN*vn_data_array[:, 2]
    Q3 = dN*vn_data_array[:, 3]
    Q4 = dN*vn_data_array[:, 4]
    Q5 = dN*vn_data_array[:, 5]
    Q6 = dN*vn_data_array[:, 6]

    # two-particle correlation
    N2_weight = dN*(dN - 1.)
    Q1_2 = abs(Q1)**2. - dN
    Q2_2 = abs(Q2)**2. - dN
    Q3_2 = abs(Q3)**2. - dN

    # four-particle correlation
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q1_4 = ((abs(Q1)**4.) - 2.*real(Q2*conj(Q1)*conj(Q1))
             - 4.*(dN - 2.)*(abs(Q1)**2.) + abs(Q2)**2.
             + 2*dN*(dN - 3.))
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))
    Q3_4 = ((abs(Q3)**4.) - 2.*real(Q6*conj(Q3)*conj(Q3))
             - 4.*(dN - 2.)*(abs(Q3)**2.) + abs(Q6)**2.
             + 2*dN*(dN - 3.))

    # C_n{4}
    C_1_4 = mean(Q1_4/N4_weight) - 2.*((mean(Q1_2/N2_weight))**2.)
    stat_err_1 = std(Q1_4/N4_weight)/sqrt(nev)
    stat_err_2 = std(Q1_2/N2_weight)/sqrt(nev)
    C_1_4_err = sqrt(stat_err_1**2.
                     + (4.*(mean(Q1_2/N2_weight))*stat_err_2)**2.)
    v1_4 = 0.0
    v1_4_err = 0.0
    if C_1_4 < 0:
        v1_4 = (-C_1_4)**0.25
        v1_4_err = 0.25*((-C_1_4)**(-0.75))*C_1_4_err

    
    C_2_4 = mean(Q2_4/N4_weight) - 2.*((mean(Q2_2/N2_weight))**2.)
    stat_err_1 = std(Q2_4/N4_weight)/sqrt(nev)
    stat_err_2 = std(Q2_2/N2_weight)/sqrt(nev)
    C_2_4_err = sqrt(stat_err_1**2.
                     + (4.*(mean(Q1_2/N2_weight))*stat_err_2)**2.)
    v2_4 = 0.0
    v2_4_err = 0.0
    if C_2_4 < 0:
        v2_4 = (-C_2_4)**0.25
        v2_4_err = 0.25*((-C_2_4)**(-0.75))*C_2_4_err

    C_3_4 = mean(Q3_4/N4_weight) - 2.*((mean(Q3_2/N2_weight))**2.)
    stat_err_1 = std(Q3_4/N4_weight)/sqrt(nev)
    stat_err_2 = std(Q3_2/N2_weight)/sqrt(nev)
    C_3_4_err = sqrt(stat_err_1**2.
                     + (4.*(mean(Q3_2/N2_weight))*stat_err_2)**2.)
    v3_4 = 0.0
    v3_4_err = 0.0
    if C_3_4 < 0:
        v3_4 = (-C_3_4)**0.25
        v3_4_err = 0.25*((-C_3_4)**(-0.75))*C_3_4_err

    results = [v1_4, v1_4_err, C_1_4, C_1_4_err,
               v2_4, v2_4_err, C_2_4, C_2_4_err,
               v3_4, v3_4_err, C_3_4, C_3_4_err,]
    return(results)


def calculate_vn4_over_vn2(vn_data_array):
    """
        this funciton computes the ratio of
        the 4-particle cumulant vn{4} over the 2-particle cumulant vn{2}
        and Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2))

            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)
            vn{2} = (<v_n*conj(v_n)>)**(1/2)

        we will use Jackknife resampling method to estimate
        the statistical error
    """
    vn_data_array = array(vn_data_array)
    nev = len(vn_data_array[:, 0])
    dN = real(vn_data_array[:, 0])
    Q1 = dN*vn_data_array[:, 1]
    Q2 = dN*vn_data_array[:, 2]
    Q3 = dN*vn_data_array[:, 3]
    Q4 = dN*vn_data_array[:, 4]
    Q5 = dN*vn_data_array[:, 5]
    Q6 = dN*vn_data_array[:, 6]

    # two-particle correlation
    N2_weight = dN*(dN - 1.)
    Q1_2 = abs(Q1)**2. - dN
    Q2_2 = abs(Q2)**2. - dN
    Q3_2 = abs(Q3)**2. - dN

    # four-particle correlation
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q1_4 = ((abs(Q1)**4.) - 2.*real(Q2*conj(Q1)*conj(Q1))
             - 4.*(dN - 2.)*(abs(Q1)**2.) + abs(Q2)**2.
             + 2*dN*(dN - 3.))
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))
    Q3_4 = ((abs(Q3)**4.) - 2.*real(Q6*conj(Q3)*conj(Q3))
             - 4.*(dN - 2.)*(abs(Q3)**2.) + abs(Q6)**2.
             + 2*dN*(dN - 3.))

    # calcualte observables with Jackknife resampling method
    r1_array = zeros(nev)
    r2_array = zeros(nev)
    r3_array = zeros(nev)
    F1_array = zeros(nev)
    F2_array = zeros(nev)
    F3_array = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        # C_n{4}
        C_1_4 = (mean(Q1_4[array_idx]/N4_weight[array_idx])
                 - 2.*((mean(Q1_2[array_idx]/N2_weight[array_idx]))**2.))
        C_1_2 = mean(Q1_2[array_idx]/N2_weight[array_idx])
        if C_1_4 < 0. and C_1_2 > 0.:
            v1_4 = (-C_1_4)**0.25
            v1_2 = sqrt(C_1_2)
            r1_array[iev] = v1_4/(v1_2 + 1e-15)
            F1_array[iev] = sqrt((v1_2**2. - v1_4**2.)
                                 /(v1_2**2. + v1_4**2. + 1e-15))

        C_2_4 = (mean(Q2_4[array_idx]/N4_weight[array_idx])
                 - 2.*((mean(Q2_2[array_idx]/N2_weight[array_idx]))**2.))
        C_2_2 = mean(Q2_2[array_idx]/N2_weight[array_idx])
        if C_2_4 < 0. and C_2_2 > 0.:
            v2_4 = (-C_2_4)**0.25
            v2_2 = sqrt(C_2_2)
            r2_array[iev] = v2_4/v2_2
            F2_array[iev] = sqrt((v2_2**2. - v2_4**2.)
                                 /(v2_2**2. + v2_4**2. + 1e-15))

        C_3_4 = (mean(Q3_4[array_idx]/N4_weight[array_idx])
                 - 2.*((mean(Q3_2[array_idx]/N2_weight[array_idx]))**2.))
        C_3_2 = mean(Q3_2[array_idx]/N2_weight[array_idx])
        if C_3_4 < 0. and C_3_2 > 0.:
            v3_4 = (-C_3_4)**0.25
            v3_2 = sqrt(C_3_2)
            r3_array[iev] = v3_4/v3_2
            F3_array[iev] = sqrt((v3_2**2. - v3_4**2.)
                                 /(v3_2**2. + v3_4**2. + 1e-15))

    r1_mean = mean(r1_array)
    r1_err = sqrt((nev - 1.)/nev*sum((r1_array - r1_mean)**2.))
    r2_mean = mean(r2_array)
    r2_err = sqrt((nev - 1.)/nev*sum((r2_array - r2_mean)**2.))
    r3_mean = mean(r3_array)
    r3_err = sqrt((nev - 1.)/nev*sum((r3_array - r3_mean)**2.))
    
    F1_mean = mean(F1_array)
    F1_err = sqrt((nev - 1.)/nev*sum((F1_array - F1_mean)**2.))
    F2_mean = mean(F2_array)
    F2_err = sqrt((nev - 1.)/nev*sum((F2_array - F2_mean)**2.))
    F3_mean = mean(F3_array)
    F3_err = sqrt((nev - 1.)/nev*sum((F3_array - F3_mean)**2.))

    results = [r1_mean, r1_err, F1_mean, F1_err,
               r2_mean, r2_err, F2_mean, F2_err,
               r3_mean, r3_err, F3_mean, F3_err]
    return(results)


def calculate_vn6_over_vn4(vn_data_array):
    """
        this funciton computes the ratio of
        the 6-particle cumulant vn{6} over the 4-particle cumulant vn{4}
            cn{6} = <<6>> - 9<<2>><<4>> + 12<<2>>^3
            vn{6} = (cn{6}/4)**(1/6)
            vn{4} = (2 <v_n*conj(v_n)>**2 - <(v_n*conj(v_n))**2.>)**(1/4)
        and compute skewness estimator gamma_1
            gamma_1 = -6\sqrt{2}*vn{4}^2*(vn{4} - vn{6})
                                         /(vn{2}^2 - vn{4}^2)^(3/2)
        we will use Jackknife resampling method to estimate
        the statistical error
    """
    vn_data_array = array(vn_data_array)
    nev = len(vn_data_array[:, 0])
    dN = real(vn_data_array[:, 0])
    Q1 = dN*vn_data_array[:, 1]
    Q2 = dN*vn_data_array[:, 2]
    Q3 = dN*vn_data_array[:, 3]
    Q4 = dN*vn_data_array[:, 4]
    Q5 = dN*vn_data_array[:, 5]
    Q6 = dN*vn_data_array[:, 6]

    # two-particle correlation
    N2_weight = dN*(dN - 1.)
    Q2_2 = abs(Q2)**2. - dN

    # four-particle correlation
    N4_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)
    Q2_4 = ((abs(Q2)**4.) - 2.*real(Q4*conj(Q2)*conj(Q2))
             - 4.*(dN - 2.)*(abs(Q2)**2.) + abs(Q4)**2.
             + 2*dN*(dN - 3.))

    # six-particle correlation
    N6_weight = dN*(dN - 1.)*(dN - 2.)*(dN - 3.)*(dN - 4.)*(dN - 5.)
    Q2_6 = (abs(Q2)**6. + 9*(abs(Q4)**2.)*(abs(Q2)**2.)
            - 6.*real(Q4*Q2*conj(Q2)*conj(Q2)*conj(Q2))
            + 4.*real(Q6*conj(Q2)*conj(Q2)*conj(Q2))
            - 12.*real(Q6*conj(Q4)*conj(Q2))
            + 18.*(dN - 4.)*real(Q4*conj(Q2)*conj(Q2))
            + 4.*(abs(Q6)**2.)
            - 9.*(dN - 4.)*((abs(Q2)**4.) + (abs(Q4)**2.))
            + 18.*(dN - 5.)*(dN - 2.)*(abs(Q2)**2.)
            - 6.*dN*(dN - 4.)*(dN - 5.))

    # calcualte observables with Jackknife resampling method
    r2_array = zeros(nev)
    gamma1_array = zeros(nev)
    for iev in range(nev):
        array_idx = [True]*nev
        array_idx[iev] = False
        array_idx = array(array_idx)

        # C_n{4}
        C_2_2 = mean(Q2_2[array_idx]/N2_weight[array_idx])
        C_2_4 = (mean(Q2_4[array_idx]/N4_weight[array_idx])
                 - 2.*(C_2_2**2.))
        C_2_6 = (mean(Q2_6[array_idx]/N6_weight[array_idx])
                 - 9.*C_2_2*mean(Q2_4[array_idx]/N4_weight[array_idx])
                 + 12.*(C_2_2**3.))
        if C_2_6 > 0. and C_2_4 < 0. and C_2_2 > 0.:
            v2_2 = sqrt(C_2_2)
            v2_6 = (C_2_6/4.)**(1./6.)
            v2_4 = (-C_2_4)**(1./4.)
            r2_array[iev] = v2_6/v2_4
            gamma1_array[iev] = (-6.*sqrt(2)*(v2_4**2.)*(v2_4 - v2_6)
                                 /(v2_2**2. - v2_4**2.)**(1.5))

    r2_mean = mean(r2_array)
    r2_err = sqrt((nev - 1.)/nev*sum((r2_array - r2_mean)**2.))
    gamma1_mean = mean(gamma1_array)
    gamma1_err = sqrt((nev - 1.)/nev*sum((gamma1_array - gamma1_mean)**2.))

    results = [r2_mean, r2_err, gamma1_mean, gamma1_err]
    return(results)


def calculate_vn_eta(dN_array, vn_array):
    nev, neta   = dN_array.shape
    dN_array    = dN_array.reshape((nev, 1, neta))
    vn_ref      = sum(dN_array*vn_array, axis=2)/(sum(dN_array, axis=2) + 1e-15)
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
        vn_SP          = mean(vn_SP_ev[array_idx, :, :], axis=0)/sqrt(vn_den)
        vn_SP_array[iev, :, :] = vn_SP
    vn_SP_mean = mean(vn_SP_array, axis=0)
    vn_SP_err  = sqrt((nev - 1.)/nev*sum((vn_SP_array - vn_SP_mean)**2., axis=0))
    return([vn_SP_mean, vn_SP_err])


def calculate_rn_eta(eta_array, dN_array, vn_array):
    nev, neta = dN_array.shape
    dN_array  = dN_array.reshape((nev, 1, neta))
    Qn_array  = vn_array
    Qnshape   = Qn_array.shape
    nQn       = Qnshape[1]
    
    # calculate the reference flow vector for every event
    eta_b_min    = 2.5
    eta_b_max    = 4.0
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
            Qn1_interp = interp(eta_ref1_tmp, eta_array, Qn_array[iev, iorder, :])
            Qn2_interp = interp(eta_ref2_tmp, eta_array, Qn_array[iev, iorder, :])
            Qn_ref1_vec.append(sum(dN1_interp*Qn1_interp)/(sum(dN1_interp) + 1e-15))
            Qn_ref2_vec.append(sum(dN2_interp*Qn2_interp)/(sum(dN2_interp) + 1e-15))
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
    return([rn_mean, rn_err, rnn_mean, rnn_err])


file_folder_list = glob(path.join(working_folder, '*'))
nev = len(file_folder_list)
for ipart, particle_id in enumerate(particle_list):
    print("processing %s ..." % particle_name_list[ipart])
    
    # first particle yield dN/dy
    if particle_id == '9999':
        file_name = 'particle_%s_vndata_eta_-0.5_0.5.dat' % particle_id
    else:
        file_name = 'particle_%s_vndata_y_-0.5_0.5.dat' % particle_id

    dN_dy = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))
        temp_data = nan_to_num(temp_data)

        dN_dy.append(temp_data[0, 1])

    dN_dy = array(dN_dy)
    dN_dy_avg = mean(dN_dy)
    dN_dy_avg_err = std(dN_dy)/sqrt(nev)

    # then <pT>, vn, dN/(2pi dy pT dpT), vn{SP}(pT)
    if particle_id == '9999':
        file_name = 'particle_%s_vndata_diff_eta_-0.5_0.5.dat' % particle_id
    else:
        file_name = 'particle_%s_vndata_diff_y_-0.5_0.5.dat' % particle_id
    file_name_ref = 'particle_9999_vndata_diff_eta_0.5_2.5.dat'
   
    pT_array = []
    dN_array = []
    vn_phenix_array = []
    vn_star_array = []
    vn_alice_array = []
    vn_cms_array = []
    vn_cms_arrays_for_rn = []
    vn_atlas_array = []
    vn_diff_phenix_real = []; vn_diff_phenix_imag = [];
    vn_diff_phenix_denorm1 = []; vn_diff_phenix_denorm2 = []
    vn_diff_star_real = []; vn_diff_star_imag = [];
    vn_diff_star_denorm1 = []; vn_diff_star_denorm2 = []
    vn_diff_alice_real = []; vn_diff_alice_imag = [];
    vn_diff_alice_denorm1 = []; vn_diff_alice_denorm2 = []
    vn_diff_2PC_real = []; vn_diff_2PC_imag = []; vn_diff_2PC_denorm = []
    vn_diff_cms_real = []; vn_diff_cms_imag = [];
    vn_diff_cms_denorm1 = []; vn_diff_cms_denorm2 = []
    vn_diff_atlas_real = []; vn_diff_atlas_imag = [];
    vn_diff_atlas_denorm1 = []; vn_diff_atlas_denorm2 = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))
        temp_data = nan_to_num(temp_data)
        temp_data_ref = loadtxt(path.join(results_folder, file_name_ref))
        temp_data_ref = nan_to_num(temp_data_ref)
        
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
        if particle_id == "9999":
            temp_vn_arrays = (
                    calculate_vn_arrays_for_rn_ratios(temp_data))
            vn_cms_arrays_for_rn.append(temp_vn_arrays)
        
        # vn with ATLAS pT cut
        temp_vn_array = calcualte_inte_vn(0.5, 3.0, temp_data)
        vn_atlas_array.append(temp_vn_array)

        # pT-differential vn using scalar-product method
        # vn{SP}(pT) with PHENIX pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff, temp_vn2 = (
                        calculate_diff_vn_single_event(0.2, 2.0, temp_data,
                                                       temp_data_ref))
        vn_diff_phenix_real.append(temp_vn_diff_real);
        vn_diff_phenix_imag.append(temp_vn_diff_imag);
        vn_diff_phenix_denorm1.append(temp_dn_diff);
        vn_diff_phenix_denorm2.append(temp_vn2);

        # vn{SP}(pT) with STAR pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff, temp_vn2 = (
                        calculate_diff_vn_single_event(0.15, 2.0, temp_data,
                                                       temp_data_ref))
        vn_diff_star_real.append(temp_vn_diff_real);
        vn_diff_star_imag.append(temp_vn_diff_imag);
        vn_diff_star_denorm1.append(temp_dn_diff);
        vn_diff_star_denorm2.append(temp_vn2);
        
        # vn{SP}(pT) with ALICE pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff, temp_vn2 = (
                        calculate_diff_vn_single_event(0.2, 3.0, temp_data,
                                                       temp_data_ref))
        vn_diff_alice_real.append(temp_vn_diff_real);
        vn_diff_alice_imag.append(temp_vn_diff_imag);
        vn_diff_alice_denorm1.append(temp_dn_diff);
        vn_diff_alice_denorm2.append(temp_vn2);
        
        # vn{SP}(pT) with CMS pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff, temp_vn2 = (
                        calculate_diff_vn_single_event(0.3, 3.0, temp_data,
                                                       temp_data_ref))
        vn_diff_cms_real.append(temp_vn_diff_real);
        vn_diff_cms_imag.append(temp_vn_diff_imag);
        vn_diff_cms_denorm1.append(temp_dn_diff);
        vn_diff_cms_denorm2.append(temp_vn2);
        
        # vn{SP}(pT) with ATLAS pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff, temp_vn2 = (
                        calculate_diff_vn_single_event(0.5, 3.0, temp_data,
                                                       temp_data_ref))
        vn_diff_atlas_real.append(temp_vn_diff_real);
        vn_diff_atlas_imag.append(temp_vn_diff_imag);
        vn_diff_atlas_denorm1.append(temp_dn_diff);
        vn_diff_atlas_denorm2.append(temp_vn2);
        
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
    pT_interp = linspace(0.15, 2.95, 30)
    dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra+1e-30)))
    dN_interp_err = interp(pT_interp, pT_spectra, dN_spectra_err)
    mean_pT_1 = sum(pT_interp**2.*dN_interp)/sum(pT_interp*dN_interp)
    mean_pT_1_upper = (sum(pT_interp**2.*(dN_interp+dN_interp_err))
                     /sum(pT_interp*(dN_interp+dN_interp_err)))
    mean_pT_1_lower = (sum(pT_interp**2.*(dN_interp-dN_interp_err))
                     /sum(pT_interp*(dN_interp-dN_interp_err)))
    mean_pT_1_err = max(abs(mean_pT_1_upper - mean_pT_1), 
                      abs(mean_pT_1 - mean_pT_1_lower))

    # calcualte vn{2}
    vn_phenix_2, vn_phenix_2_err = calcualte_vn_2(vn_phenix_array)
    vn_star_2, vn_star_2_err = calcualte_vn_2(vn_star_array)
    vn_alice_2, vn_alice_2_err = calcualte_vn_2(vn_alice_array)
    vn_cms_2, vn_cms_2_err = calcualte_vn_2(vn_cms_array)
    vn_atlas_2, vn_atlas_2_err = calcualte_vn_2(vn_atlas_array)

    if (particle_id == '9999'):
        vn_alice_array2 = array(vn_alice_array)
        vn_cms_array2 = array(vn_cms_array)
        vn_atlas_array2 = array(vn_atlas_array)
        # calculate non-linear response coefficents with ALICE pT cut
        nonlinear_response_alice = calculate_nonlinear_reponse(vn_alice_array2)
        # calculate non-linear response coefficents with CMS pT cut
        nonlinear_response_cms = calculate_nonlinear_reponse(vn_cms_array2)
        # calculate non-linear response coefficents with ATLAS pT cut
        nonlinear_response_atlas = calculate_nonlinear_reponse(vn_atlas_array2)
        
        # calculate symmetric cumulant coefficents with ALICE pT cut
        SC_alice = calculate_symmetric_cumulant(vn_alice_array)

        # calculate vn{4}
        vn4_alice = calculate_vn4(vn_alice_array)
        vn4_cms = calculate_vn4(vn_cms_array)
        vn4_atlas = calculate_vn4(vn_atlas_array)
        # calculate vn{4}/vn{2} and vn{6}/vn{4}
        vn4_over_vn2_alice = calculate_vn4_over_vn2(vn_alice_array)
        vn4_over_vn2_cms = calculate_vn4_over_vn2(vn_cms_array)
        vn4_over_vn2_atlas = calculate_vn4_over_vn2(vn_atlas_array)
        vn6_over_vn4_alice = calculate_vn6_over_vn4(vn_alice_array)
        vn6_over_vn4_cms = calculate_vn6_over_vn4(vn_cms_array)
        vn6_over_vn4_atlas = calculate_vn6_over_vn4(vn_atlas_array)

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
            vn_diff_phenix_real, vn_diff_phenix_imag,
            vn_diff_phenix_denorm1, vn_diff_phenix_denorm2)

    vn_diff_SP_star, vn_diff_SP_star_err = calculate_vn_diff_SP(
            vn_diff_star_real, vn_diff_star_imag,
            vn_diff_star_denorm1, vn_diff_star_denorm2)
    
    vn_diff_SP_alice, vn_diff_SP_alice_err = calculate_vn_diff_SP(
            vn_diff_alice_real, vn_diff_alice_imag,
            vn_diff_alice_denorm1, vn_diff_alice_denorm2)
    
    vn_diff_SP_cms, vn_diff_SP_cms_err = calculate_vn_diff_SP(
            vn_diff_cms_real, vn_diff_cms_imag,
            vn_diff_cms_denorm1, vn_diff_cms_denorm2)
        
    vn_diff_SP_atlas, vn_diff_SP_atlas_err = calculate_vn_diff_SP(
            vn_diff_atlas_real, vn_diff_atlas_imag,
            vn_diff_atlas_denorm1, vn_diff_atlas_denorm2)
    
    # calcualte vn[2](pT)
    vn_diff_2PC, vn_diff_2PC_err = calculate_vn_diff_2PC(
            vn_diff_2PC_real, vn_diff_2PC_imag, vn_diff_2PC_denorm)
    
    # then particle rapidity distribution
    if particle_id == '9999':
        file_name = 'particle_%s_dNdeta_pT_0.2_3.dat' % particle_id
    else:
        file_name = 'particle_%s_dNdy_pT_0.2_3.dat' % particle_id

    eta_array = []
    dN_array = []
    vn_array = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))
        temp_data = nan_to_num(temp_data)

        eta_array.append(temp_data[:, 0])
        dN_array.append(temp_data[:, 1])
        temp_vn_array = []
        for iorder in range(1, n_order):
            vn_real = temp_data[:, 6*iorder-3]
            vn_imag = temp_data[:, 6*iorder-1]
            vn = vn_real + 1j*vn_imag
            temp_vn_array.append(vn)
        vn_array.append(temp_vn_array)

    eta_array = array(eta_array)
    dN_array = array(dN_array)
    vn_array = array(vn_array)

    eta_point = mean(eta_array, 0)
    dNdeta = mean(dN_array, 0)
    dNdeta_err = std(dN_array, 0)/sqrt(nev)
    vn_SP_eta, vn_SP_eta_err = calculate_vn_eta(dN_array, vn_array)
    rn_eta, rn_eta_err, rnn_eta, rnn_eta_err = calculate_rn_eta(eta_point, dN_array, vn_array)
    vn_eta_real = mean(real(vn_array), 0)
    vn_eta_real_err = std(real(vn_array), 0)/sqrt(nev)
   
    ###########################################################################
    # finally, output all the results
    ###########################################################################
    
    if (particle_id =='9999'):
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
        
        # output symmetric cumulants for ALICE pt cut
        output_filename = ("symmetric_cumulant_ALICE.dat")
        f = open(output_filename, 'w')
        f.write("# type  value  stat. err\n")
        for i in range(len(symmetric_cumulant_name_list)):
            f.write("%s  %.10e  %.10e\n"
                    % (symmetric_cumulant_name_list[i],
                       SC_alice[2*i], SC_alice[2*i+1]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output vn4 for ALICE pt cut
        output_filename = ("charged_hadron_vn4_ALICE.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{4}  vn{4}_err  Cn{4}  Cn{4}_err\n")
        for i in range(1, 4):
            f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                    % (i, vn4_alice[4*i-4], vn4_alice[4*i-3],
                       vn4_alice[4*i-2], vn4_alice[4*i-1]))
        f.close()
        shutil.move(output_filename, avg_folder)

        # output vn4 for CMS pt cut
        output_filename = ("charged_hadron_vn4_CMS.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{4}  vn{4}_err  Cn{4}  Cn{4}_err\n")
        for i in range(1, 4):
            f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                    % (i, vn4_cms[4*i-4], vn4_cms[4*i-3],
                       vn4_cms[4*i-2], vn4_cms[4*i-1]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output vn4 for ATLAS pt cut
        output_filename = ("charged_hadron_vn4_ATLAS.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{4}  vn{4}_err  Cn{4}  Cn{4}_err\n")
        for i in range(1, 4):
            f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                    % (i, vn4_atlas[4*i-4], vn4_atlas[4*i-3],
                       vn4_atlas[4*i-2], vn4_atlas[4*i-1]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output vn4/vn2 ratio for ALICE pt cut
        output_filename = ("charged_hadron_vn4_over_vn2_ALICE.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{4}/vn{2}  (vn{4}/vn{2})_err  Fn  Fn_err \n")
        f.write("# Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2)) \n")
        for i in range(1, 4):
            f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                    % (i, vn4_over_vn2_alice[4*i-4], vn4_over_vn2_alice[4*i-3],
                       vn4_over_vn2_alice[4*i-2], vn4_over_vn2_alice[4*i-1]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output vn4/vn2 ratio for CMS pt cut
        output_filename = ("charged_hadron_vn4_over_vn2_CMS.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{4}/vn{2}  (vn{4}/vn{2})_err  Fn  Fn_err \n")
        f.write("# Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2)) \n")
        for i in range(1, 4):
            f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                    % (i, vn4_over_vn2_cms[4*i-4], vn4_over_vn2_cms[4*i-3],
                       vn4_over_vn2_cms[4*i-2], vn4_over_vn2_cms[4*i-1]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output vn4/vn2 ratio for ATLAS pt cut
        output_filename = ("charged_hadron_vn4_over_vn2_ATLAS.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{4}/vn{2}  (vn{4}/vn{2})_err  Fn  Fn_err \n")
        f.write("# Fn = sqrt((vn{2}^2 - vn{4}^2)/(vn{2}^2 + vn{4}^2)) \n")
        for i in range(1, 4):
            f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                    % (i, vn4_over_vn2_atlas[4*i-4], vn4_over_vn2_atlas[4*i-3],
                       vn4_over_vn2_atlas[4*i-2], vn4_over_vn2_atlas[4*i-1]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output vn6/vn4 ratio for ALICE pt cut
        output_filename = ("charged_hadron_vn6_over_vn4_ALICE.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{6}/vn{4}  (vn{6}/vn{4})_err  gamma_1  gamma_1_err \n")
        f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                % (2, vn6_over_vn4_alice[0], vn6_over_vn4_alice[1],
                   vn6_over_vn4_alice[2], vn6_over_vn4_alice[3]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output vn6/vn4 ratio for CMS pt cut
        output_filename = ("charged_hadron_vn6_over_vn4_CMS.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{6}/vn{4}  (vn{6}/vn{4})_err  gamma_1  gamma_1_err \n")
        f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                % (2, vn6_over_vn4_cms[0], vn6_over_vn4_cms[1],
                   vn6_over_vn4_cms[2], vn6_over_vn4_cms[3]))
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output vn6/vn4 ratio for ATLAS pt cut
        output_filename = ("charged_hadron_vn6_over_vn4_ATLAS.dat")
        f = open(output_filename, 'w')
        f.write("# n  vn{6}/vn{4}  (vn{6}/vn{4})_err  gamma_1  gamma_1_err \n")
        f.write("%d  %.10e  %.10e  %.10e  %.10e\n"
                % (2, vn6_over_vn4_atlas[0], vn6_over_vn4_atlas[1],
                   vn6_over_vn4_atlas[2], vn6_over_vn4_atlas[3]))
        f.close()
        shutil.move(output_filename, avg_folder)

    output_filename = ("%s_integrated_observables.dat"
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    f.write("dN/dy= %.10e +/- %.10e\n" % (dN_dy_avg, dN_dy_avg_err))
    f.write("<pT>= %.10e +/- %.10e\n" % (mean_pT, mean_pT_err))
    f.write("<pT(>0.15)>= %.10e +/- %.10e\n" % (mean_pT_1, mean_pT_1_err))
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
                       % particle_name_list[ipart])
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
                       % particle_name_list[ipart])
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
                       % particle_name_list[ipart])
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
                       % particle_name_list[ipart])
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
                       % particle_name_list[ipart])
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
                       % particle_name_list[ipart])
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

    output_filename = ("%s_rapidity_distribution.dat" 
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    if(particle_id == '9999'):
        f.write("#eta  dN/deta  dN/deta_err  vn{2}(eta)  vn{2}(eta)_err\n")
    else:
        f.write("#y  dN/dy  dN/dy_err  vn{2}(y)  vn{2}(y)_err\n")
    for ieta in range(len(eta_point)):
        f.write("%.10e  %.10e  %.10e  "
                % (eta_point[ieta], dNdeta[ieta], dNdeta_err[ieta]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  %.10e  %.10e  "
                    % (vn_SP_eta[iorder-1, ieta],
                       vn_SP_eta_err[iorder-1, ieta],
                       vn_eta_real[iorder-1, ieta],
                       vn_eta_real_err[iorder-1, ieta]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    
    if (particle_id == '9999'):
        output_filename = ("%s_vn_distribution_PHENIX.dat"
                           % particle_name_list[ipart])
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
                           % particle_name_list[ipart])
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
                           % particle_name_list[ipart])
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
                           % particle_name_list[ipart])
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
                           % particle_name_list[ipart])
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
                           % (particle_name_list[ipart],
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
                                       % (particle_name_list[ipart],
                                          pT_trig[ipTtrig],
                                          pT_trig[ipTtrig+1]))
                    f = open(output_filename, 'w')
                    f.write("#pT_mid  rn  rn_err (n = 2, 3, 4)\n")

        output_filename = ("%s_rn_eta.dat" 
                           % particle_name_list[ipart])
        f = open(output_filename, 'w')
        f.write("#eta  rn(eta)  rn_err(eta)  rnn(eta)  rnn_err(eta)\n")
        for ieta in range(len(eta_point)-1):
            f.write("%.10e  " % eta_point[ieta])
            for iorder in range(0, n_order-1):
                f.write("%.10e  %.10e  %.10e  %.10e  "
                        % (rn_eta[iorder, ieta],
                           rn_eta_err[iorder, ieta],
                           rnn_eta[iorder, ieta],
                           rnn_eta_err[iorder, ieta]))
            f.write("\n")
        f.close()
        shutil.move(output_filename, avg_folder)
        
        # output flow event-plane correlation
        output_filename = ("%s_event_plane_correlation_ALICE.dat"
                           % particle_name_list[ipart])
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
                           % particle_name_list[ipart])
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

print("Analysis is done.")

