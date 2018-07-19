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

rap_region = "-0.5_0.5"

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


def calcualte_vn_2(vn_data_array):
    """
        this function computes vn{2} and its stat. err.
        self correlation is substracted
    """
    vn_data_array = array(vn_data_array)
    nev = len(vn_data_array[:, 0])
    dN = vn_data_array[:, 0]
    dN = dN.reshape(len(dN), 1)
    Qn_array = dN*vn_data_array[:, 1:]
    corr = 1./(dN*(dN - 1.))*(Qn_array*conj(Qn_array) - dN)
    vn_2 = sqrt(real(mean(corr, 0))) + 1e-30
    vn_2_err = std(real(corr), 0)/sqrt(nev)/2./vn_2
    return(vn_2, vn_2_err)


def calcualte_vn_2_delta_eta(vn_A_data_array, vn_B_data_array):
    """
        this function computes vn{2} and its stat. err.
        sample A and B are assumed to have no overlap
        no self correlation substraction needed
    """
    vn_A_data_array = array(vn_A_data_array)
    vn_B_data_array = array(vn_B_data_array)
    nev = len(vn_A_data_array[:, 0])
    dN_A = vn_A_data_array[:, 0]
    dN_A = dN_A.reshape(len(dN_A), 1)
    dN_B = vn_B_data_array[:, 0]
    dN_B = dN_B.reshape(len(dN_B), 1)
    Qn_A_array = dN_A*vn_A_data_array[:, 1:]
    Qn_B_array = dN_B*vn_B_data_array[:, 1:]
    corr = 1./(dN_A*dN_B)*Qn_A_array*conj(Qn_B_array)
    vn_2 = sqrt(real(mean(corr, 0))) + 1e-30
    vn_2_err = std(real(corr), 0)/sqrt(nev)/2./vn_2
    return(vn_2, vn_2_err)

def calcualte_three_plane_correlations(
    vn_data_array_A, vn_data_array_B, vn_data_array_ref, type):
    """
        this function compute three plane correlation
        vn_data_array is a matrix [event_idx, vn_order]
        type == 0:
            we assume no overlap among samples A, B, and ref
        type == 1:
            for opposite sign particle pairs (A, B)
            samples A and B do not have any overlap in this case
            we assume samples A, B \elem sample ref
        type == 2:
            for same sign particle pairs (A, B)
            samples A = sample B in this case
            we assume samples A, B \elem sample ref
    """
    vn_data_array_A = array(vn_data_array_A)
    vn_data_array_B = array(vn_data_array_B)
    vn_data_array_ref = array(vn_data_array_ref)
    nev = len(vn_data_array_ref[:, 0])
    dN_A = vn_data_array_A[:, 0]
    dN_B = vn_data_array_B[:, 0]
    dN_ref = vn_data_array_ref[:, 0]
    dN_A = dN_A.reshape(len(dN_A), 1)
    dN_B = dN_B.reshape(len(dN_B), 1)
    dN_ref = dN_A.reshape(len(dN_ref), 1)
    Qn_array_A = dN_A*vn_data_array_A[:, 1:]
    Qn_array_B = dN_B*vn_data_array_B[:, 1:]
    Qn_array_ref = dN_ref*vn_data_array_ref[:, 1:]
    Q1_A = Qn_array_A[:, 0]
    Q2_A = Qn_array_A[:, 1]
    Q3_A = Qn_array_A[:, 2]
    Q4_A = Qn_array_A[:, 3]
    Q5_A = Qn_array_A[:, 4]
    Q6_A = Qn_array_A[:, 5]
    Q1_B = Qn_array_B[:, 0]
    Q2_B = Qn_array_B[:, 1]
    Q3_B = Qn_array_B[:, 2]
    Q4_B = Qn_array_B[:, 3]
    Q5_B = Qn_array_B[:, 4]
    Q6_B = Qn_array_B[:, 5]
    Q1_ref = Qn_array_ref[:, 0]
    Q2_ref = Qn_array_ref[:, 1]
    Q3_ref = Qn_array_ref[:, 2]
    Q4_ref = Qn_array_ref[:, 3]
    Q5_ref = Qn_array_ref[:, 4]
    Q6_ref = Qn_array_ref[:, 5]

    # C_112 = <Re{(V_1)(V_1)conj(V_2)}>
    if type == 0:
        corr_112_data = 1./(dN_A*dN_B*dN_ref)*Q1_A*Q1_B*conj(Q2_ref)
    elif type == 1:
        corr_112_data = 1./(dN_A*dN_B*(dN_ref - 2.))*(
            Q1_A*Q1_B*conj(Q2_ref) - Q1_A*conj(Q1_B) - Q1_B*conj(Q1_A)
        )
    elif type == 2:
        corr_112_data = 1./(dN_A*(dN_B - 1.)*(dN_ref - 2.))*(
            Q1_A*Q1_B*conj(Q2_ref) - Q1_A*conj(Q1_B) - Q1_B*conj(Q1_A)
            - Q2_A*conj(Q2_ref) + 2.*dN_A
        )
    corr_112 = mean(real(corr_112_data))
    corr_112_err = std(real(corr_112_data))/sqrt(nev)
    
    # C_123 = <Re{(V_1)(V_2)conj(V_3)}>
    if type == 0:
        corr_123_data = 1./(dN_A*dN_B*dN_ref)*Q1_A*Q2_B*conj(Q3_ref)
    elif type == 1:
        corr_123_data = 1./(dN_A*dN_B*(dN_ref - 2.))*(
            Q1_A*Q2_B*conj(Q3_ref) - Q1_A*conj(Q1_B) - Q2_B*conj(Q2_A)
        )
    elif type == 2:
        corr_123_data = 1./(dN_A*(dN_B - 1.)*(dN_ref - 2.))*(
            Q1_A*Q2_B*conj(Q3_ref) - Q1_A*conj(Q1_B) - Q2_B*conj(Q2_A)
            - Q3_A*conj(Q3_ref) + 2.*dN_A
        )
    corr_123 = mean(real(corr_123_data))
    corr_123_err = std(real(corr_123_data))/sqrt(nev)

    # C_224 = <Re{(V_2)(V_2)conj(V_4)}>
    if type == 0:
        corr_224_data = 1./(dN_A*dN_B*dN_ref)*Q2_A*Q2_B*conj(Q4_ref)
    elif type == 1:
        corr_224_data = 1./(dN_A*dN_B*(dN_ref - 2.))*(
            Q2_A*Q2_B*conj(Q4_ref) - Q2_A*conj(Q2_B) - Q2_B*conj(Q2_A)
        )
    elif type == 2:
        corr_224_data = 1./(dN_A*(dN_B - 1.)*(dN_ref - 2.))*(
            Q2_A*Q2_B*conj(Q4_ref) - Q2_A*conj(Q2_B) - Q2_B*conj(Q2_A)
            - Q4_A*conj(Q4_ref) + 2.*dN_A
        )
    corr_224 = mean(real(corr_224_data))
    corr_224_err = std(real(corr_224_data))/sqrt(nev)

    # C_235 = <Re{(V_2)(V_3)conj(V_5)}>
    if type == 0:
        corr_235_data = 1./(dN_A*dN_B*dN_ref)*Q2_A*Q3_B*conj(Q5_ref)
    elif type == 1:
        corr_235_data = 1./(dN_A*dN_B*(dN_ref - 2.))*(
            Q2_A*Q3_B*conj(Q5_ref) - Q2_A*conj(Q2_B) - Q3_B*conj(Q3_A)
        )
    elif type == 2:
        corr_235_data = 1./(dN_A*(dN_B - 1.)*(dN_ref - 2.))*(
            Q2_A*Q3_B*conj(Q5_ref) - Q2_A*conj(Q2_B) - Q3_B*conj(Q3_A)
            - Q5_A*conj(Q5_ref) + 2.*dN_A
        )
    corr_235 = mean(real(corr_235_data))
    corr_235_err = std(real(corr_235_data))/sqrt(nev)
    
    # C_134 = <Re{(V_1)(V_3)conj(V_4)}>
    if type == 0:
        corr_134_data = 1./(dN_A*dN_B*dN_ref)*Q1_A*Q3_B*conj(Q4_ref)
    elif type == 1:
        corr_134_data = 1./(dN_A*dN_B*(dN_ref - 2.))*(
            Q1_A*Q3_B*conj(Q4_ref) - Q1_A*conj(Q1_B) - Q3_B*conj(Q3_A)
        )
    elif type == 2:
        corr_134_data = 1./(dN_A*(dN_B - 1.)*(dN_ref - 2.))*(
            Q1_A*Q3_B*conj(Q4_ref) - Q1_A*conj(Q1_B) - Q3_B*conj(Q3_A)
            - Q4_A*conj(Q4_ref) + 2.*dN_A
        )
    corr_134 = mean(real(corr_134_data))
    corr_134_err = std(real(corr_134_data))/sqrt(nev)

    # C_246 = <Re{(V_2)(V_4)conj(V_6)}>
    if type == 0:
        corr_246_data = 1./(dN_A*dN_B*dN_ref)*Q2_A*Q4_B*conj(Q6_ref)
    elif type == 1:
        corr_246_data = 1./(dN_A*dN_B*(dN_ref - 2.))*(
            Q2_A*Q4_B*conj(Q6_ref) - Q2_A*conj(Q2_B) - Q4_B*conj(Q4_A)
        )
    elif type == 2:
        corr_246_data = 1./(dN_A*(dN_B - 1.)*(dN_ref - 2.))*(
            Q2_A*Q4_B*conj(Q6_ref) - Q2_A*conj(Q2_B) - Q4_B*conj(Q4_A)
            - Q6_A*conj(Q6_ref) + 2.*dN_A
        )
    corr_246 = mean(real(corr_246_data))
    corr_246_err = std(real(corr_246_data))/sqrt(nev)

    # C_336 = <Re{(V_3)(V_3)conj(V_6)}>
    if type == 0:
        corr_336_data = 1./(dN_A*dN_B*dN_ref)*Q3_A*Q3_B*conj(Q6_ref)
    elif type == 1:
        corr_336_data = 1./(dN_A*dN_B*(dN_ref - 2.))*(
            Q3_A*Q3_B*conj(Q6_ref) - Q3_A*conj(Q3_B) - Q3_B*conj(Q3_A)
        )
    elif type == 2:
        corr_336_data = 1./(dN_A*(dN_B - 1.)*(dN_ref - 2.))*(
            Q3_A*Q3_B*conj(Q6_ref) - Q3_A*conj(Q3_B) - Q3_B*conj(Q3_A)
            - Q6_A*conj(Q6_ref) + 2.*dN_A
        )
    corr_336 = mean(real(corr_336_data))
    corr_336_err = std(real(corr_336_data))/sqrt(nev)

    results = [corr_112, corr_123, corr_224, corr_235, corr_134, corr_246, corr_336]
    results_err = [corr_112_err, corr_123_err, corr_224_err, corr_235_err, corr_134_err, corr_246_err, corr_336_err]
    return(results, results_err)


file_folder_list = glob(path.join(working_folder, '*'))
nev = len(file_folder_list)
print("processing three particle correlations ...")
# load the file
file_name_pos = 'particle_9998_vndata_diff_eta_%s.dat' % rap_region
file_name_neg = 'particle_-9998_vndata_diff_eta_%s.dat' % rap_region
file_name_ch = 'particle_9999_vndata_diff_eta_%s.dat' % rap_region
file_name_ch_1 = 'particle_9999_vndata_diff_eta_-2_-0.5.dat'
file_name_ch_2 = 'particle_9999_vndata_diff_eta_0.5_2.dat'
   
vn_pos_phenix_array = []; vn_neg_phenix_array = []; vn_ch_phenix_array = []
vn_pos_star_array = []; vn_neg_star_array = []; vn_ch_star_array = []
vn_pos_alice_array = []; vn_neg_alice_array = []; vn_ch_alice_array = []
vn_pos_cms_array = []; vn_neg_cms_array = []; vn_ch_cms_array = []
vn_pos_atlas_array = []; vn_neg_atlas_array = []; vn_ch_atlas_array = []
vn_ch_cms_array_1 = []; vn_ch_cms_array_2 = []
for ifolder in range(nev):
    results_folder = path.abspath(file_folder_list[ifolder])
    temp_data_ch = loadtxt(path.join(results_folder, file_name_ch))
    temp_data_pos = loadtxt(path.join(results_folder, file_name_pos))
    temp_data_neg = loadtxt(path.join(results_folder, file_name_neg))
    temp_data_ch_1 = loadtxt(path.join(results_folder, file_name_ch_1))
    temp_data_ch_2 = loadtxt(path.join(results_folder, file_name_ch_2))
        
    # pT-integrated vn
    # vn with PHENIX pT cut
    temp_vn_array = calcualte_inte_vn(0.2, 2.0, temp_data_ch)
    vn_ch_phenix_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.2, 2.0, temp_data_pos)
    vn_pos_phenix_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.2, 2.0, temp_data_neg)
    vn_neg_phenix_array.append(temp_vn_array)

    # vn with STAR pT cut
    temp_vn_array = calcualte_inte_vn(0.2, 2.0, temp_data_ch)
    vn_ch_star_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.2, 2.0, temp_data_pos)
    vn_pos_star_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.2, 2.0, temp_data_neg)
    vn_neg_star_array.append(temp_vn_array)

    # vn with ALICE pT cut
    temp_vn_array = calcualte_inte_vn(0.2, 3.0, temp_data_ch)
    vn_ch_alice_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.2, 3.0, temp_data_pos)
    vn_pos_alice_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.2, 3.0, temp_data_neg)
    vn_neg_alice_array.append(temp_vn_array)
        
    # vn with CMS pT cut
    temp_vn_array = calcualte_inte_vn(0.3, 3.0, temp_data_ch)
    vn_ch_cms_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.3, 3.0, temp_data_pos)
    vn_pos_cms_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.3, 3.0, temp_data_neg)
    vn_neg_cms_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.3, 3.0, temp_data_ch_1)
    vn_ch_cms_array_1.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.3, 3.0, temp_data_ch_2)
    vn_ch_cms_array_2.append(temp_vn_array)
        
    # vn with ATLAS pT cut
    temp_vn_array = calcualte_inte_vn(0.5, 3.0, temp_data_ch)
    vn_ch_atlas_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.5, 3.0, temp_data_pos)
    vn_pos_atlas_array.append(temp_vn_array)
    temp_vn_array = calcualte_inte_vn(0.5, 3.0, temp_data_neg)
    vn_neg_atlas_array.append(temp_vn_array)

# now we perform event average
# calcualte vn{2}
vn_phenix_2, vn_phenix_2_err = calcualte_vn_2(vn_ch_phenix_array)
vn_star_2, vn_star_2_err = calcualte_vn_2(vn_ch_star_array)
vn_alice_2, vn_alice_2_err = calcualte_vn_2(vn_ch_alice_array)
vn_cms_2, vn_cms_2_err = calcualte_vn_2(vn_ch_cms_array)
#print(vn_cms_2, vn_cms_2_err)
#vn_cms_2, vn_cms_2_err = calcualte_vn_2_delta_eta(vn_ch_cms_array_1, vn_ch_cms_array_2)
#print(vn_cms_2, vn_cms_2_err)
vn_atlas_2, vn_atlas_2_err = calcualte_vn_2(vn_ch_atlas_array)

# calculate flow event-plane correlation
corr_os_cms, corr_os_cms_err = calcualte_three_plane_correlations(
            vn_pos_cms_array, vn_neg_cms_array, vn_ch_cms_array, 1)
corr_ss_cms, corr_ss_cms_err = calcualte_three_plane_correlations(
            vn_pos_cms_array, vn_pos_cms_array, vn_ch_cms_array, 2)
corr_os_star, corr_os_star_err = calcualte_three_plane_correlations(
            vn_pos_star_array, vn_neg_star_array, vn_ch_star_array, 1)
corr_pp_star, corr_pp_star_err = calcualte_three_plane_correlations(
            vn_pos_star_array, vn_pos_star_array, vn_ch_star_array, 2)
corr_mm_star, corr_mm_star_err = calcualte_three_plane_correlations(
            vn_neg_star_array, vn_neg_star_array, vn_ch_star_array, 2)
corr_ch_star, corr_ch_star_err = calcualte_three_plane_correlations(
            vn_ch_star_array, vn_ch_star_array, vn_ch_star_array, 2)

###########################################################################
# finally, output all the results
###########################################################################
    
# output charged hadron flow coefficients
output_filename = ("charged_hadron_vn_STAR.dat")
f = open(output_filename, 'w')
f.write("#v1{2}  v1{2}_err  v2{2}  v2{2}_err  v3{2}  v3{2}_err  "
        +"v4{2}  v4{2}_err  v5{2}  v5{2}_err  v6{2}  v6{2}_err  "
        +"v7{2}  v7{2}_err\n")
for i in range(n_order-1):
    f.write("%.5e  %.5e  " % (vn_star_2[i], vn_star_2_err[i]))
f.close()
shutil.move(output_filename, avg_folder)

# output flow event-plane correlation
output_filename = ("three_plane_correlation_os_STAR.dat")
f = open(output_filename, 'w')
f.write("#v2{2}  v2{2}_err  C_112  C_112_err  C_123  C_123_err  "
        +"C_224  C_224_err  C_235  C_235_err  C_134  C_134_err  "
        +"C_246  C_246_err  C_336  C_336_err\n")
f.write("%.5e  %.5e  " % (vn_star_2[1], vn_star_2_err[1]))
for i in range(len(corr_os_star)):
    f.write("%.5e  %.5e  " % (corr_os_star[i], corr_os_star_err[i]))
f.close()
shutil.move(output_filename, avg_folder)

output_filename = ("three_plane_correlation_pp_STAR.dat")
f = open(output_filename, 'w')
f.write("#v2{2}  v2{2}_err  C_112  C_112_err  C_123  C_123_err  "
        +"C_224  C_224_err  C_235  C_235_err  C_134  C_134_err  "
        +"C_246  C_246_err  C_336  C_336_err\n")
f.write("%.5e  %.5e  " % (vn_star_2[1], vn_star_2_err[1]))
for i in range(len(corr_pp_star)):
    f.write("%.5e  %.5e  " % (corr_pp_star[i], corr_pp_star_err[i]))
f.close()
shutil.move(output_filename, avg_folder)

output_filename = ("three_plane_correlation_mm_STAR.dat")
f = open(output_filename, 'w')
f.write("#v2{2}  v2{2}_err  C_112  C_112_err  C_123  C_123_err  "
        +"C_224  C_224_err  C_235  C_235_err  C_134  C_134_err  "
        +"C_246  C_246_err  C_336  C_336_err\n")
f.write("%.5e  %.5e  " % (vn_star_2[1], vn_star_2_err[1]))
for i in range(len(corr_mm_star)):
    f.write("%.5e  %.5e  " % (corr_mm_star[i], corr_mm_star_err[i]))
f.close()
shutil.move(output_filename, avg_folder)

output_filename = ("three_plane_correlation_ch_STAR.dat")
f = open(output_filename, 'w')
f.write("#v2{2}  v2{2}_err  C_112  C_112_err  C_123  C_123_err  "
        +"C_224  C_224_err  C_235  C_235_err  C_134  C_134_err  "
        +"C_246  C_246_err  C_336  C_336_err\n")
f.write("%.5e  %.5e  " % (vn_star_2[1], vn_star_2_err[1]))
for i in range(len(corr_ch_star)):
    f.write("%.5e  %.5e  " % (corr_ch_star[i], corr_ch_star_err[i]))
f.close()
shutil.move(output_filename, avg_folder)

output_filename = ("three_plane_correlation_os_CMS.dat")
f = open(output_filename, 'w')
f.write("#v2{2}  v2{2}_err  C_112  C_112_err  C_123  C_123_err  "
        +"C_224  C_224_err  C_235  C_235_err  C_134  C_134_err  "
        +"C_246  C_246_err  C_336  C_336_err\n")
f.write("%.5e  %.5e  " % (vn_cms_2[1], vn_cms_2_err[1]))
for i in range(len(corr_os_cms)):
    f.write("%.5e  %.5e  " % (corr_os_cms[i], corr_os_cms_err[i]))
f.close()
shutil.move(output_filename, avg_folder)

output_filename = ("three_plane_correlation_pp_CMS.dat")
f = open(output_filename, 'w')
f.write("#v2{2}  v2{2}_err  C_112  C_112_err  C_123  C_123_err  "
        +"C_224  C_224_err  C_235  C_235_err  C_134  C_134_err  "
        +"C_246  C_246_err  C_336  C_336_err\n")
f.write("%.5e  %.5e  " % (vn_cms_2[1], vn_cms_2_err[1]))
for i in range(len(corr_ss_cms)):
    f.write("%.5e  %.5e  " % (corr_ss_cms[i], corr_ss_cms_err[i]))
f.close()
shutil.move(output_filename, avg_folder)

print("Analysis is done.")

