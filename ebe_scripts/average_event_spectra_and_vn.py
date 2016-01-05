#! /usr/bin/env python
"""
     This script performs event averaging for particle spectra and vn
     calculated from event-by-event simulations
"""

from sys import argv, exit
from os import path
from glob import glob
from numpy import *

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

try:
    working_folder = path.abspath(argv[1])
except(IndexError):
    print("Usage: average_event_spectra_and_vn.py working_folder")
    exit(1)

file_folder_list = glob(path.join(working_folder, '*'))
particle_lists = ['9999', '211', '-211', '321', '-321', '2212', '-2212']

weight_type = 0     # 0: weight with particle multiplicity, 1: weight with 1

nev = len(file_folder_list)
for ipart in range(len(particle_lists)):
    # pT-integrated quantities
    file_name = 'particle_%s_vndata.dat' % particle_lists[ipart]
    dn_array = []
    v1_real_array = []; v1_imag_array = []
    v2_real_array = []; v2_imag_array = []
    v3_real_array = []; v3_imag_array = []
    v4_real_array = []; v4_imag_array = []
    v5_real_array = []; v5_imag_array = []
    
    # read in data
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print "read in data event %d : %s ..." % (ifolder, file_name)
        temp_data = loadtxt(path.join(results_folder, file_name))

        # dN/dy
        dn_array.append(temp_data[0, 1])
        # vn
        v1_real_array.append(temp_data[1, 1])
        v1_imag_array.append(temp_data[1, 3])
        v2_real_array.append(temp_data[2, 1])
        v2_imag_array.append(temp_data[2, 3])
        v3_real_array.append(temp_data[3, 1])
        v3_imag_array.append(temp_data[3, 3])
        v4_real_array.append(temp_data[4, 1])
        v4_imag_array.append(temp_data[4, 3])
        v5_real_array.append(temp_data[5, 1])
        v5_imag_array.append(temp_data[5, 3])

    dn_array = array(dn_array)
    v1_real_array = array(v1_real_array) 
    v1_imag_array = array(v1_imag_array)
    v2_real_array = array(v2_real_array)
    v2_imag_array = array(v2_imag_array)
    v3_real_array = array(v3_real_array)
    v3_imag_array = array(v3_imag_array)
    v4_real_array = array(v4_real_array)
    v4_imag_array = array(v4_imag_array)
    v5_real_array = array(v5_real_array)
    v5_imag_array = array(v5_imag_array)

    dn_avg = mean(dn_array)
    dn_avg_err = std(dn_array)/sqrt(nev)

    if weight_type == 0:
        weight = dn_array
    elif weight_type == 1:
        weight = ones(len(dn_array))

    vn_mag_sq_array = v1_real_array**2. + v1_imag_array**2.
    v1_avg = mean(sqrt(vn_mag_sq_array)*weight)/mean(weight)
    v1_SP = sqrt(mean(vn_mag_sq_array*weight*weight)/mean(weight**2.))
    v1_avg_err = sqrt(
        (std(sqrt(vn_mag_sq_array)*weight)/mean(weight)/sqrt(nev))**2.
      + (mean(sqrt(vn_mag_sq_array)*weight)*std(weight)/(mean(weight)**2.)
        /sqrt(nev))**2.
    )
    v1_SP_err = sqrt(
        (
            std((vn_mag_sq_array)*weight*weight)
            /(2.*sqrt(mean((vn_mag_sq_array)*weight*weight)*mean(weight**2.)))
        )**2. 
      + (
            sqrt(mean(vn_mag_sq_array*weight**2.))*std(weight**2.)
            /(mean(weight**2.)**1.5)
        )**2.
    )
    
    vn_mag_sq_array = v2_real_array**2. + v2_imag_array**2.
    v2_avg = mean(sqrt(vn_mag_sq_array)*weight)/mean(weight)
    v2_SP = sqrt(mean(vn_mag_sq_array*weight*weight)/mean(weight**2.))
    v2_avg_err = sqrt(
        (std(sqrt(vn_mag_sq_array)*weight)/mean(weight)/sqrt(nev))**2.
      + (mean(sqrt(vn_mag_sq_array)*weight)*std(weight)/(mean(weight)**2.)
        /sqrt(nev))**2.
    )
    v2_SP_err = sqrt(
        (
            std((vn_mag_sq_array)*weight*weight)
            /(2.*sqrt(mean((vn_mag_sq_array)*weight*weight)*mean(weight**2.)))
        )**2. 
      + (
            sqrt(mean(vn_mag_sq_array*weight**2.))*std(weight**2.)
            /(mean(weight**2.)**1.5)
        )**2.
    )
    
    vn_mag_sq_array = v3_real_array**2. + v3_imag_array**2.
    v3_avg = mean(sqrt(vn_mag_sq_array)*weight)/mean(weight)
    v3_SP = sqrt(mean(vn_mag_sq_array*weight*weight)/mean(weight**2.))
    v3_avg_err = sqrt(
        (std(sqrt(vn_mag_sq_array)*weight)/mean(weight)/sqrt(nev))**2.
      + (mean(sqrt(vn_mag_sq_array)*weight)*std(weight)/(mean(weight)**2.)
        /sqrt(nev))**2.
    )
    v3_SP_err = sqrt(
        (
            std((vn_mag_sq_array)*weight*weight)
            /(2.*sqrt(mean((vn_mag_sq_array)*weight*weight)*mean(weight**2.)))
        )**2. 
      + (
            sqrt(mean(vn_mag_sq_array*weight**2.))*std(weight**2.)
            /(mean(weight**2.)**1.5)
        )**2.
    )

    vn_mag_sq_array = v4_real_array**2. + v4_imag_array**2.
    v4_avg = mean(sqrt(vn_mag_sq_array)*weight)/mean(weight)
    v4_SP = sqrt(mean(vn_mag_sq_array*weight*weight)/mean(weight**2.))
    v4_avg_err = sqrt(
        (std(sqrt(vn_mag_sq_array)*weight)/mean(weight)/sqrt(nev))**2.
      + (mean(sqrt(vn_mag_sq_array)*weight)*std(weight)/(mean(weight)**2.)
        /sqrt(nev))**2.
    )
    v4_SP_err = sqrt(
        (
            std((vn_mag_sq_array)*weight*weight)
            /(2.*sqrt(mean((vn_mag_sq_array)*weight*weight)*mean(weight**2.)))
        )**2. 
      + (
            sqrt(mean(vn_mag_sq_array*weight**2.))*std(weight**2.)
            /(mean(weight**2.)**1.5)
        )**2.
    )

    vn_mag_sq_array = v5_real_array**2. + v5_imag_array**2.
    v5_avg = mean(sqrt(vn_mag_sq_array)*weight)/mean(weight)
    v5_SP = sqrt(mean(vn_mag_sq_array*weight*weight)/mean(weight**2.))
    v5_avg_err = sqrt(
        (std(sqrt(vn_mag_sq_array)*weight)/mean(weight)/sqrt(nev))**2.
      + (mean(sqrt(vn_mag_sq_array)*weight)*std(weight)/(mean(weight)**2.)
        /sqrt(nev))**2.
    )
    v5_SP_err = sqrt(
        (
            std((vn_mag_sq_array)*weight*weight)
            /(2.*sqrt(mean((vn_mag_sq_array)*weight*weight)*mean(weight**2.)))
        )**2. 
      + (
            sqrt(mean(vn_mag_sq_array*weight**2.))*std(weight**2.)
            /(mean(weight**2.)**1.5)
        )**2.
    )
    
    output = []
    output.append([0, dn_avg, dn_avg_err, 0, 0])
    output.append([1, v1_avg, v1_avg_err, v1_SP, v1_SP_err])
    output.append([2, v2_avg, v2_avg_err, v2_SP, v2_SP_err])
    output.append([3, v3_avg, v3_avg_err, v3_SP, v3_SP_err])
    output.append([4, v4_avg, v4_avg_err, v4_SP, v4_SP_err])
    output.append([5, v5_avg, v5_avg_err, v5_SP, v5_SP_err])

    savetxt('pT_integrated_dN_vn_%s.dat' % particle_lists[ipart], 
            output, fmt='%.10e', delimiter='  ',)

    # pT-differential flow
    file_name_ref = 'particle_9999_vndata.dat'
    dn_ref_array = []
    v1_real_ref_array = []; v1_imag_ref_array = []
    v2_real_ref_array = []; v2_imag_ref_array = []
    v3_real_ref_array = []; v3_imag_ref_array = []
    v4_real_ref_array = []; v4_imag_ref_array = []
    v5_real_ref_array = []; v5_imag_ref_array = []

    file_name = 'particle_%s_vndata_diff.dat' % particle_lists[ipart]
    pT_array = []
    dn_pT_array = []
    v1_real_pT_array = []; v1_imag_pT_array = []
    v2_real_pT_array = []; v2_imag_pT_array = []
    v3_real_pT_array = []; v3_imag_pT_array = []
    v4_real_pT_array = []; v4_imag_pT_array = []
    v5_real_pT_array = []; v5_imag_pT_array = []
    
    # read in data
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print "read in data event %d : %s ..." % (ifolder, file_name)
        temp_data_ref = loadtxt(path.join(results_folder, file_name_ref))
        temp_data = loadtxt(path.join(results_folder, file_name))

        # dN/dy
        dn_ref_array.append(temp_data_ref[0, 1])
        # vn
        v1_real_ref_array.append(temp_data_ref[1, 1])
        v1_imag_ref_array.append(temp_data_ref[1, 3])
        v2_real_ref_array.append(temp_data_ref[2, 1])
        v2_imag_ref_array.append(temp_data_ref[2, 3])
        v3_real_ref_array.append(temp_data_ref[3, 1])
        v3_imag_ref_array.append(temp_data_ref[3, 3])
        v4_real_ref_array.append(temp_data_ref[4, 1])
        v4_imag_ref_array.append(temp_data_ref[4, 3])
        v5_real_ref_array.append(temp_data_ref[5, 1])
        v5_imag_ref_array.append(temp_data_ref[5, 3])
        
        pT_array = temp_data[:, 0]
        # dN/(2pi dy pT dpT)
        dn_pT_array.append(temp_data[:, 2])
        # vn
        v1_real_pT_array.append(temp_data[:, 4])
        v1_imag_pT_array.append(temp_data[:, 6])
        v2_real_pT_array.append(temp_data[:, 8])
        v2_imag_pT_array.append(temp_data[:, 10])
        v3_real_pT_array.append(temp_data[:, 12])
        v3_imag_pT_array.append(temp_data[:, 14])
        v4_real_pT_array.append(temp_data[:, 16])
        v4_imag_pT_array.append(temp_data[:, 18])
        v5_real_pT_array.append(temp_data[:, 20])
        v5_imag_pT_array.append(temp_data[:, 22])
    
    dn_ref_array = array(dn_ref_array)
    v1_real_ref_array = array(v1_real_ref_array) 
    v1_imag_ref_array = array(v1_imag_ref_array)
    v2_real_ref_array = array(v2_real_ref_array)
    v2_imag_ref_array = array(v2_imag_ref_array)
    v3_real_ref_array = array(v3_real_ref_array)
    v3_imag_ref_array = array(v3_imag_ref_array)
    v4_real_ref_array = array(v4_real_ref_array)
    v4_imag_ref_array = array(v4_imag_ref_array)
    v5_real_ref_array = array(v5_real_ref_array)
    v5_imag_ref_array = array(v5_imag_ref_array)
    v1_real_ref_array = v1_real_ref_array.reshape(len(v1_real_ref_array), 1)
    v1_imag_ref_array = v1_imag_ref_array.reshape(len(v1_imag_ref_array), 1)
    v2_real_ref_array = v2_real_ref_array.reshape(len(v2_real_ref_array), 1)
    v2_imag_ref_array = v2_imag_ref_array.reshape(len(v2_imag_ref_array), 1)
    v3_real_ref_array = v3_real_ref_array.reshape(len(v3_real_ref_array), 1)
    v3_imag_ref_array = v3_imag_ref_array.reshape(len(v3_imag_ref_array), 1)
    v4_real_ref_array = v4_real_ref_array.reshape(len(v4_real_ref_array), 1)
    v4_imag_ref_array = v4_imag_ref_array.reshape(len(v4_imag_ref_array), 1)
    v5_real_ref_array = v5_real_ref_array.reshape(len(v5_real_ref_array), 1)
    v5_imag_ref_array = v5_imag_ref_array.reshape(len(v5_imag_ref_array), 1)
    
    dn_pT_array = array(dn_pT_array)
    v1_real_pT_array = array(v1_real_pT_array) 
    v1_imag_pT_array = array(v1_imag_pT_array)
    v2_real_pT_array = array(v2_real_pT_array)
    v2_imag_pT_array = array(v2_imag_pT_array)
    v3_real_pT_array = array(v3_real_pT_array)
    v3_imag_pT_array = array(v3_imag_pT_array)
    v4_real_pT_array = array(v4_real_pT_array)
    v4_imag_pT_array = array(v4_imag_pT_array)
    v5_real_pT_array = array(v5_real_pT_array)
    v5_imag_pT_array = array(v5_imag_pT_array)
    
    dn_pT_avg = mean(dn_pT_array, 0)
    dn_pT_avg_err = std(dn_pT_array, 0)/sqrt(nev)

    if weight_type == 0:
        weight_pT = dn_pT_array
        weight_ref = dn_ref_array
        weight_ref = weight_ref.reshape(len(weight_ref), 1)
    elif weight_type == 1:
        weight_pT = ones(dn_pT_array.shape)
        weight_ref = ones(len(dn_ref_array), 1)
    
    vn_pT_mag_sq_array = v1_real_pT_array**2. + v1_imag_pT_array**2.
    v1_avg_pT = mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)/mean(weight_pT, 0)
    numerator = (
        (v1_real_pT_array*weight_pT)*(v1_real_ref_array*weight_ref)
        - (v1_imag_pT_array*weight_pT)*(v1_imag_ref_array*weight_ref))
    denominator = (v1_real_ref_array**2. + v1_imag_ref_array**2.)*weight_ref**2.
    v1_SP_pT = mean(numerator, 0)/mean(sqrt(denominator))
    v1_avg_pT_err = sqrt(
        (std(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)
         /mean(weight_pT, 0)/sqrt(nev))**2.
        + (mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)*std(weight_pT, 0)
           /(mean(weight_pT, 0)**2.))**2.
    )
    v1_SP_pT_err = sqrt(
        (std(numerator, 0)/mean(sqrt(denominator))/sqrt(nev))**2.
        + (mean(numerator, 0)*std(sqrt(denominator))
          /(mean(sqrt(denominator)))**2.)**2.
    )
    
    vn_pT_mag_sq_array = v2_real_pT_array**2. + v2_imag_pT_array**2.
    v2_avg_pT = mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)/mean(weight_pT, 0)
    numerator = (
        (v2_real_pT_array*weight_pT)*(v2_real_ref_array*weight_ref)
        - (v2_imag_pT_array*weight_pT)*(v2_imag_ref_array*weight_ref))
    denominator = (v2_real_ref_array**2. + v2_imag_ref_array**2.)*weight_ref**2.
    v2_SP_pT = mean(numerator, 0)/mean(sqrt(denominator))
    v2_avg_pT_err = sqrt(
        (std(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)
         /mean(weight_pT, 0)/sqrt(nev))**2.
        + (mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)*std(weight_pT, 0)
           /(mean(weight_pT, 0)**2.))**2.
    )
    v2_SP_pT_err = sqrt(
        (std(numerator, 0)/mean(sqrt(denominator))/sqrt(nev))**2.
        + (mean(numerator, 0)*std(sqrt(denominator))
          /(mean(sqrt(denominator)))**2.)**2.
    )
    
    vn_pT_mag_sq_array = v3_real_pT_array**2. + v3_imag_pT_array**2.
    v3_avg_pT = mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)/mean(weight_pT, 0)
    numerator = (
        (v3_real_pT_array*weight_pT)*(v3_real_ref_array*weight_ref)
        - (v3_imag_pT_array*weight_pT)*(v3_imag_ref_array*weight_ref))
    denominator = (v3_real_ref_array**2. + v3_imag_ref_array**2.)*weight_ref**2.
    v3_SP_pT = mean(numerator, 0)/mean(sqrt(denominator))
    v3_avg_pT_err = sqrt(
        (std(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)
         /mean(weight_pT, 0)/sqrt(nev))**2.
        + (mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)*std(weight_pT, 0)
           /(mean(weight_pT, 0)**2.))**2.
    )
    v3_SP_pT_err = sqrt(
        (std(numerator, 0)/mean(sqrt(denominator))/sqrt(nev))**2.
        + (mean(numerator, 0)*std(sqrt(denominator))
          /(mean(sqrt(denominator)))**2.)**2.
    )
    
    vn_pT_mag_sq_array = v4_real_pT_array**2. + v4_imag_pT_array**2.
    v4_avg_pT = mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)/mean(weight_pT, 0)
    numerator = (
        (v4_real_pT_array*weight_pT)*(v4_real_ref_array*weight_ref)
        - (v4_imag_pT_array*weight_pT)*(v4_imag_ref_array*weight_ref))
    denominator = (v4_real_ref_array**2. + v4_imag_ref_array**2.)*weight_ref**2.
    v4_SP_pT = mean(numerator, 0)/mean(sqrt(denominator))
    v4_avg_pT_err = sqrt(
        (std(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)
         /mean(weight_pT, 0)/sqrt(nev))**2.
        + (mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)*std(weight_pT, 0)
           /(mean(weight_pT, 0)**2.))**2.
    )
    v4_SP_pT_err = sqrt(
        (std(numerator, 0)/mean(sqrt(denominator))/sqrt(nev))**2.
        + (mean(numerator, 0)*std(sqrt(denominator))
          /(mean(sqrt(denominator)))**2.)**2.
    )
    
    vn_pT_mag_sq_array = v5_real_pT_array**2. + v5_imag_pT_array**2.
    v5_avg_pT = mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)/mean(weight_pT, 0)
    numerator = (
        (v5_real_pT_array*weight_pT)*(v5_real_ref_array*weight_ref)
        - (v5_imag_pT_array*weight_pT)*(v5_imag_ref_array*weight_ref))
    denominator = (v5_real_ref_array**2. + v5_imag_ref_array**2.)*weight_ref**2.
    v5_SP_pT = mean(numerator, 0)/mean(sqrt(denominator))
    v5_avg_pT_err = sqrt(
        (std(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)
         /mean(weight_pT, 0)/sqrt(nev))**2.
        + (mean(sqrt(vn_pT_mag_sq_array)*weight_pT, 0)*std(weight_pT, 0)
           /(mean(weight_pT, 0)**2.))**2.
    )
    v5_SP_pT_err = sqrt(
        (std(numerator, 0)/mean(sqrt(denominator))/sqrt(nev))**2.
        + (mean(numerator, 0)*std(sqrt(denominator))
          /(mean(sqrt(denominator)))**2.)**2.
    )

    output_pT_avg = array([pT_array, dn_pT_avg, dn_pT_avg_err, 
                                     v1_avg_pT, v1_avg_pT_err, 
                                     v2_avg_pT, v2_avg_pT_err, 
                                     v3_avg_pT, v3_avg_pT_err, 
                                     v4_avg_pT, v4_avg_pT_err, 
                                     v5_avg_pT, v5_avg_pT_err,])
    
    output_pT_SP = array([pT_array, dn_pT_avg, dn_pT_avg_err, 
                                    v1_SP_pT, v1_SP_pT_err,
                                    v2_SP_pT, v2_SP_pT_err,
                                    v3_SP_pT, v3_SP_pT_err,
                                    v4_SP_pT, v4_SP_pT_err,
                                    v5_SP_pT, v5_SP_pT_err])

    savetxt('pT_differential_dN_vn_avg_%s.dat' % particle_lists[ipart], 
            output_pT_avg.transpose(), fmt='%.10e', delimiter='  ',)
    savetxt('pT_differential_dN_vn_SP_%s.dat' % particle_lists[ipart], 
            output_pT_SP.transpose(), fmt='%.10e', delimiter='  ',)


print "Analysis is done."
