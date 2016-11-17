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

particle_list = ['211', '-211', '321', '-321', '2212', '-2212', 
                 '3122', '-3122', '3312', '-3312', '3334', '-3334',
                 '333', '9999']
particle_name_list = ['pion_p', 'pion_m', 'kaon_p', 'kaon_m', 'proton', 
                      'anti_proton', 'Lambda', 'anti_Lambda', 'Xi_m',
                      'anti_Xi_p', 'Omega', 'anti_Omega', 'phi',
                      'charged_hadron']


file_folder_list = glob(path.join(working_folder, '*'))
nev = len(file_folder_list)
for ipart, particle_id in enumerate(particle_list):
    print("processing %s ..." % particle_name_list[ipart])
    
    # first event-by-event particle yield dN/dy distribution
    if particle_id == '9999':
        file_name = 'particle_%s_yield_distribution_eta.dat' % particle_id
    else:
        file_name = 'particle_%s_yield_distribution_y.dat' % particle_id

    N_array = []
    pN_array = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))
        
        N_event = temp_data[:, 0]
        pN_event = temp_data[:, 1]

        # record particle spectra
        if ifolder == 0:
            N_array.append(N_event)
        pN_array.append(pN_event)

    # now we perform event average
    N_array = array(N_array)
    pN_array = array(pN_array)
    pN_avg = mean(pN_array, 0)

    ###########################################################################
    # finally, output all the results
    ###########################################################################
    
    output_filename = ("%s_yield_distribution.dat" 
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    f.write("# N  P(N)\n")
    for ipT in range(len(pN_avg)):
        if (pN_avg[ipT] > 1e-10):
            f.write("%d  %.10e\n"
                     % (N_array[0, ipT], pN_avg[ipT]))
    f.close()
    shutil.move(output_filename, avg_folder)

print "Analysis is done."

