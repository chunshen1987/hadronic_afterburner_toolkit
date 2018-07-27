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
    print("Usage: %s working_folder results_folder" % argv[0])
    exit(1)

rap_region = "-1_1"

n_order = 7

file_folder_list = glob(path.join(working_folder, '*'))
nev = len(file_folder_list)

print("processing two particle correlations ...")
file_name_ch = 'particle_9999_vn2_eta_%s.dat' % rap_region
Qn2_array = []
for ifolder in range(nev):
    results_folder = path.abspath(file_folder_list[ifolder])
    temp_data_vn2 = loadtxt(path.join(results_folder, file_name_ch))
    Qn2_array.append(temp_data_vn2)
Qn2_array = array(Qn2_array)
output_filename = ("two_particle_correlation_STAR.dat")
f = open(output_filename, 'w')
f.write("# n  vn{2}^2  vn{2}^2_err\n")
Npair = mean(Qn2_array[:, 0, 3])
Npair_err = sqrt(mean(Qn2_array[:, 0, 4]) - Npair**2.)/sqrt(nev)
f.write("%s  %.5e  %.5e\n" % (0, Npair, Npair_err))
for ii in range(1, 9):
    vn2_ch = mean(Qn2_array[:, ii, 3])
    vn2_ch_err = sqrt(mean(Qn2_array[:, ii, 4]) - vn2_ch**2.)/sqrt(nev)
    vn2_ch = vn2_ch/Npair
    vn2_ch_err = vn2_ch_err/Npair
    f.write("%s  %.5e  %.5e\n" % (ii, vn2_ch, vn2_ch_err))
f.close()
shutil.move(output_filename, avg_folder)

print("processing two particle correlations delta eta dependence ...")
file_name_ch = 'particle_9999_vn2_eta12_pT_0.2_3.dat'
Qn2_array = []
for ifolder in range(nev):
    results_folder = path.abspath(file_folder_list[ifolder])
    temp_data_vn2 = loadtxt(path.join(results_folder, file_name_ch))
    Qn2_array.append(temp_data_vn2)
Qn2_array = array(Qn2_array)
output_filename = ("two_particle_correlation_delta_eta12_STAR.dat")
f = open(output_filename, 'w')
f.write("# rap  vn{2}^2  vn{2}^2_err\n")
Npair = mean(Qn2_array[:, :, 3], axis=0)
#Npair_err = sqrt(mean(Qn2_array[:, :, 4], axis=0) - Npair**2.)/sqrt(nev)
Npair_err = std(Qn2_array[:, :, 3], axis=0)/sqrt(nev)
output = []
output.append(Qn2_array[0, :, 0])
output.append(Npair)
output.append(Npair_err)
for ii in range(1, 9):
    #vn2_ch = mean(Qn2_array[:, :, 4*ii+3], axis=0)
    #vn2_ch_err = sqrt(mean(Qn2_array[:, :, 4*ii+4], axis=0)
    #                  - vn2_ch**2.)/sqrt(nev)
    #vn2_ch = vn2_ch/Npair
    #vn2_ch_err = vn2_ch_err/Npair
    vn2_ch = mean(Qn2_array[:, :, 4*ii+1]*Qn2_array[:, :, 3], axis=0)/Npair
    vn2_ch_err = sqrt(mean(Qn2_array[:, :, 4*ii+2]**2., axis=0))/sqrt(nev)
    output.append(vn2_ch)
    output.append(vn2_ch_err)
output = array(output)
output = output.transpose()
for irap in range(len(Npair)):
    f.write("%.5e  "*19 % tuple(output[irap, :]))
    f.write("\n")
f.close()
shutil.move(output_filename, avg_folder)

print("processing three particle correlations ...")
# load the file
file_name_ss = 'particle_9999_Cmnk_ss_eta_%s.dat' % rap_region
file_name_os = 'particle_9999_Cmnk_os_eta_%s.dat' % rap_region
file_name_ch = 'particle_9999_Cmnk_eta_%s.dat' % rap_region
file_name_vn2 = 'particle_9999_vn2_eta_%s.dat' % rap_region
file_name_spvn = 'particle_9999_vndata_diff_eta_%s.dat' % rap_region

C_mnk_ch_array = []
C_mnk_ss_array = []
C_mnk_os_array = []
Qn2_array = []
dN_array = []; pT_array = []
for ifolder in range(nev):
    results_folder = path.abspath(file_folder_list[ifolder])
    temp_data_ch = loadtxt(path.join(results_folder, file_name_ch))
    temp_data_ss = loadtxt(path.join(results_folder, file_name_ss))
    temp_data_os = loadtxt(path.join(results_folder, file_name_os))
    temp_data_vn2 = loadtxt(path.join(results_folder, file_name_vn2))
    temp_data_spvn = loadtxt(path.join(results_folder, file_name_spvn))

    C_mnk_ch_array.append(temp_data_ch)
    C_mnk_ss_array.append(temp_data_ss)
    C_mnk_os_array.append(temp_data_os)
    Qn2_array.append(temp_data_vn2[:, 3:5])
    dN_array.append(temp_data_spvn[:, 2])
    pT_array.append(temp_data_spvn[:, 0])

C_mnk_ch_array = array(C_mnk_ch_array)
C_mnk_ss_array = array(C_mnk_ss_array)
C_mnk_os_array = array(C_mnk_os_array)
Qn2_array = array(Qn2_array)
dN_array = array(dN_array)
pT_array = array(pT_array)
n_pT = len(pT_array[0, :])
pT_spectra = zeros([n_pT])
for ipT in range(len(pT_array[0, :])):
    dN_temp = sum(dN_array[:, ipT]*pT_array[:, ipT])
    if (dN_temp > 0):
        pT_spectra[ipT] = (
                sum(pT_array[:, ipT]**2.*dN_array[:, ipT])/dN_temp)
    else:
        pT_spectra[ipT] = mean(pT_array[:, ipT])
dN_spectra = mean(pT_array*dN_array, 0)/pT_spectra   # dN/(2pi dy pT dpT)
pT_interp = linspace(0.05, 3.5, 40)
dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra+1e-30)))
mean_pT_sq = sum(pT_interp**3.*dN_interp)/sum(pT_interp*dN_interp)
pT_interp = linspace(0.2, 3.5, 40)
dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra+1e-30)))
mean_pT = sum(pT_interp**2.*dN_interp)/sum(pT_interp*dN_interp)
mean_pT_sq_2 = sum(pT_interp**3.*dN_interp)/sum(pT_interp*dN_interp)
factor = mean_pT**2./mean_pT_sq
factor2 = mean_pT_sq_2/mean_pT_sq
print(factor, factor2)

nev = len(C_mnk_ch_array[:, 0, 1])

corr_label = ['000', '112', '123', '224', '235', '134', '246', '336', '347']
output_filename = ("three_plane_correlation_STAR.dat")
f = open(output_filename, 'w')
f.write("#C_mnk_ss  C_mnk_ss_err  C_mnk_os  C_mnk_os_err  C_mnk_ch  C_mnk_ch_err\n")
for ii in range(1, len(corr_label)):
    momentum_conservation = 0.0
    #if ii == 1:
    #    momentum_conservation = - factor*2.*Qn2_array[:, 2, 0]
    C_mnk_ch_avg = (
            sum(C_mnk_ch_array[:, 0, 1]*C_mnk_ch_array[:, ii, 1]
                + momentum_conservation)/(sum(C_mnk_ch_array[:, 0, 1])))
    C_mnk_ch_err = (
            sum(C_mnk_ch_array[:, 0, 1]*C_mnk_ch_array[:, ii, 2])
            /(sum(C_mnk_ch_array[:, 0, 1]))/sqrt(nev))
    C_mnk_os_avg = (sum(C_mnk_os_array[:, 0, 1]*C_mnk_os_array[:, ii, 1])
            /(sum(C_mnk_os_array[:, 0, 1])))
    C_mnk_os_err = (sum(C_mnk_os_array[:, 0, 1]*C_mnk_os_array[:, ii, 2])
            /(sum(C_mnk_os_array[:, 0, 1]))/sqrt(nev))
    C_mnk_ss_avg = (sum(C_mnk_ss_array[:, 0, 1]*C_mnk_ss_array[:, ii, 1])
            /(sum(C_mnk_ss_array[:, 0, 1])))
    C_mnk_ss_err = (sum(C_mnk_ss_array[:, 0, 1]*C_mnk_ss_array[:, ii, 2])
            /(sum(C_mnk_ss_array[:, 0, 1]))/sqrt(nev))
    f.write("%s  %.5e  %.5e  %.5e  %.5e  %.5e  %.5e\n"
            % (corr_label[ii], C_mnk_ss_avg, C_mnk_ss_err,
               C_mnk_os_avg, C_mnk_os_err, C_mnk_ch_avg, C_mnk_ch_err))
f.close()
shutil.move(output_filename, avg_folder)

print("processing three particle correlations as a function of delta eta ...")
eta_type_list = ['eta12', 'eta13']
for itype in range(len(eta_type_list)):
    eta_type = eta_type_list[itype]
    # load the file
    file_name_ss = 'particle_9999_Cmnk_%s_ss_pT_0.2_3.dat' % eta_type
    file_name_os = 'particle_9999_Cmnk_%s_os_pT_0.2_3.dat' % eta_type
    file_name_ch = 'particle_9999_Cmnk_%s_pT_0.2_3.dat' % eta_type
    
    C_mnk_ch_array = []
    C_mnk_ss_array = []
    C_mnk_os_array = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data_ch = loadtxt(path.join(results_folder, file_name_ch))
        temp_data_ss = loadtxt(path.join(results_folder, file_name_ss))
        temp_data_os = loadtxt(path.join(results_folder, file_name_os))
    
        C_mnk_ch_array.append(temp_data_ch)
        C_mnk_ss_array.append(temp_data_ss)
        C_mnk_os_array.append(temp_data_os)
    
    C_mnk_ch_array = array(C_mnk_ch_array)
    C_mnk_ss_array = array(C_mnk_ss_array)
    C_mnk_os_array = array(C_mnk_os_array)
    nev, nrap, n_corr = C_mnk_ch_array.shape
    
    results_ch = zeros([nrap, n_corr])
    results_ch[:, 0] = mean(C_mnk_ch_array[:, :, 0], axis=0)
    results_ss = zeros([nrap, n_corr])
    results_ss[:, 0] = mean(C_mnk_ss_array[:, :, 0], axis=0)
    results_os = zeros([nrap, n_corr])
    results_os[:, 0] = mean(C_mnk_os_array[:, :, 0], axis=0)
    for ii in range(0, len(corr_label)):
        if ii == 0:
            results_ch[:, 2*ii+1] = mean(C_mnk_ch_array[:, :, 2*ii+1], axis=0)
            results_ch[:, 2*ii+2] = (
                sqrt(mean(C_mnk_ch_array[:, :, 2*ii+2]**2., axis=0))/sqrt(nev))
            results_ss[:, 2*ii+1] = mean(C_mnk_ss_array[:, :, 2*ii+1], axis=0)
            results_ss[:, 2*ii+2] = (
                sqrt(mean(C_mnk_ss_array[:, :, 2*ii+2]**2., axis=0))/sqrt(nev))
            results_os[:, 2*ii+1] = mean(C_mnk_os_array[:, :, 2*ii+1], axis=0)
            results_os[:, 2*ii+2] = (
                sqrt(mean(C_mnk_os_array[:, :, 2*ii+2]**2., axis=0))/sqrt(nev))
        else:
            results_ch[:, 2*ii+1] = (
                mean(C_mnk_ch_array[:, :, 2*ii+1]*C_mnk_ch_array[:, :, 1], axis=0)
                /mean(C_mnk_ch_array[:, :, 1], axis=0))
            results_ch[:, 2*ii+2] = (
                sqrt(mean(C_mnk_ch_array[:, :, 2*ii+2]**2., axis=0))/sqrt(nev))
            results_ss[:, 2*ii+1] = (
                mean(C_mnk_ss_array[:, :, 2*ii+1]*C_mnk_ss_array[:, :, 1], axis=0)
                /mean(C_mnk_ss_array[:, :, 1], axis=0))
            results_ss[:, 2*ii+2] = (
                sqrt(mean(C_mnk_ss_array[:, :, 2*ii+2]**2., axis=0))/sqrt(nev))
            results_os[:, 2*ii+1] = (
                mean(C_mnk_os_array[:, :, 2*ii+1]*C_mnk_os_array[:, :, 1], axis=0)
                /mean(C_mnk_os_array[:, :, 1], axis=0))
            results_os[:, 2*ii+2] = (
                sqrt(mean(C_mnk_os_array[:, :, 2*ii+2]**2., axis=0))/sqrt(nev))
    ncol = 2*len(corr_label) + 1
    output_filename = ("three_plane_correlation_ch_delta_%s_STAR.dat" % eta_type)
    f = open(output_filename, 'w')
    f.write("# %s  C_nmk  C_nmk_err (000, 112, 123, 224, 235, 134, 246, 336, 347)\n"
            % eta_type)
    for irap in range(len(Npair)):
        f.write("%.5e  "*ncol % tuple(results_ch[irap, :]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    output_filename = ("three_plane_correlation_ss_delta_%s_STAR.dat" % eta_type)
    f = open(output_filename, 'w')
    f.write("# %s  C_nmk  C_nmk_err (000, 112, 123, 224, 235, 134, 246, 336, 347)\n"
            % eta_type)
    for irap in range(len(Npair)):
        f.write("%.5e  "*ncol % tuple(results_ss[irap, :]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    output_filename = ("three_plane_correlation_os_delta_%s_STAR.dat" % eta_type)
    f = open(output_filename, 'w')
    f.write("# %s  C_nmk  C_nmk_err (000, 112, 123, 224, 235, 134, 246, 336, 347)\n"
            % eta_type)
    for irap in range(len(Npair)):
        f.write("%.5e  "*ncol % tuple(results_os[irap, :]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)

print("processing four particle correlations ...")
# load the file
file_name_ch = 'particle_9999_SCmn_eta_%s.dat' % rap_region
SC_mn_ch_array = []
for ifolder in range(nev):
    results_folder = path.abspath(file_folder_list[ifolder])
    temp_data_ch = loadtxt(path.join(results_folder, file_name_ch))
    SC_mn_ch_array.append(temp_data_ch)
SC_mn_ch_array = array(SC_mn_ch_array)
nev = len(SC_mn_ch_array[:, 0, 1])
corr_label = ['00', '32', '42', '52', '43', '53']
output_filename = ("symmetric_cumulant_STAR.dat")
f = open(output_filename, 'w')
f.write("# name SC_mn_ch  SC_mn_ch_err\n")
for ii in range(1, len(corr_label)):
    SC_mn_ch_avg = (sum(SC_mn_ch_array[:, 0, 1]*SC_mn_ch_array[:, ii, 1])
            /(sum(SC_mn_ch_array[:, 0, 1])))
    SC_mn_ch_err = (sum(SC_mn_ch_array[:, 0, 1]*SC_mn_ch_array[:, ii, 2])
            /(sum(SC_mn_ch_array[:, 0, 1]))/sqrt(nev))
    f.write("%s  %.5e  %.5e\n" % (corr_label[ii], SC_mn_ch_avg, SC_mn_ch_err))
f.close()
shutil.move(output_filename, avg_folder)

print("Analysis is done.")

