#! /usr/bin/env python3
"""
     This script performs event averaging for the balance correlation function
     calculated from event-by-event simulations

     Format from balance correlation function file:
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
    print("Usage: %s working_folder results_folder" % sys.argv[0])
    exit(1)

file_folder_list = glob(path.join(working_folder, '*'))
pair_list = ["9998_-9998"]

nev = len(file_folder_list)
for ipart in pair_list:
    file_name = 'Balance_function_%s_Delta_y.dat' % ipart
    event_avg_data = loadtxt(
            path.join(file_folder_list[0],
                      'Balance_function_9998_-9998_Delta_y.dat'))*0.0
    event_avg_data = zeros([len(event_avg_data[:, 0]), 15])
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print("processing %s/%s ..." % (results_folder, file_name))
        temp_data = loadtxt(path.join(results_folder, file_name))
        event_avg_data[:, 0] += temp_data[:, 0]
        event_avg_data[:, 5] += temp_data[:, 3]       # rho_2(OS)
        event_avg_data[:, 6] += temp_data[:, 3]**2.
        event_avg_data[:, 7] += temp_data[:, 4]       # rho_1^2(OS)
        event_avg_data[:, 8] += temp_data[:, 4]**2.
        event_avg_data[:, 11] += temp_data[:, 6]       # rho_2(SS)
        event_avg_data[:, 12] += temp_data[:, 6]**2.
        event_avg_data[:, 13] += temp_data[:, 7]       # rho_1^2(SS)
        event_avg_data[:, 14] += temp_data[:, 7]**2.

    event_avg_data = event_avg_data/nev
    # rho_2(OS)_err
    event_avg_data[:, 6] = sqrt(event_avg_data[:, 6]
                                - event_avg_data[:, 5]**2.)/sqrt(nev)
    # rho_1^2(OS)_err
    event_avg_data[:, 8] = sqrt(event_avg_data[:, 8]
                                - event_avg_data[:, 7]**2.)/sqrt(nev)
    # rho_2(SS)_err
    event_avg_data[:, 12] = sqrt(event_avg_data[:, 12]
                                - event_avg_data[:, 11]**2.)/sqrt(nev)
    # rho_1^2(SS)_err
    event_avg_data[:, 14] = sqrt(event_avg_data[:, 14]
                                - event_avg_data[:, 13]**2.)/sqrt(nev)
    # C_2(OS)
    norm = sum(event_avg_data[:, 7])/sum(event_avg_data[:, 5])
    event_avg_data[:, 3] = norm*event_avg_data[:, 5]/event_avg_data[:, 7]
    event_avg_data[:, 4] = norm*sqrt(
          (event_avg_data[:, 6]/event_avg_data[:, 7])**2.
        + (event_avg_data[:, 5]*event_avg_data[:, 8]/event_avg_data[:, 7]**2.)**2.)
    # C_2(SS)
    norm = sum(event_avg_data[:, 13])/sum(event_avg_data[:, 11])
    event_avg_data[:, 9] = norm*event_avg_data[:, 11]/event_avg_data[:, 13]
    event_avg_data[:, 10] = norm*sqrt(
          (event_avg_data[:, 12]/event_avg_data[:, 13])**2.
        + (event_avg_data[:, 11]*event_avg_data[:, 14]/event_avg_data[:, 13]**2.)**2.)
    # Delta C_2
    event_avg_data[:, 1] = event_avg_data[:, 3] - event_avg_data[:, 9]
    event_avg_data[:, 2] = sqrt(event_avg_data[:, 4]**2.
                                + event_avg_data[:, 10]**2.)

    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ',
            header=("DeltaY  Delta C_2  Delta C_2_err  "
                    + "C_2(OS)  C_2(OS)_err  rho_2(OS)  rho_2(OS)_err  "
                    + "rho_1^2(OS)  rho_1^2(OS)_err  "
                    + "C_2(SS)  C_2(SS)_err  rho_2(SS)  rho_2(SS)_err  "
                    + "rho_1^2(SS)  rho_1^2(SS)_err"))
    shutil.move(file_name, avg_folder)
    
    file_name = 'Balance_function_%s_Delta_phi.dat' % ipart
    event_avg_data = loadtxt(
            path.join(file_folder_list[0],
                      'Balance_function_9998_-9998_Delta_phi.dat'))*0.0
    event_avg_data = zeros([len(event_avg_data[:, 0]), 15])
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print("processing %s/%s ..." % (results_folder, file_name))
        temp_data = loadtxt(path.join(results_folder, file_name))
        event_avg_data[:, 0] += temp_data[:, 0]
        event_avg_data[:, 5] += temp_data[:, 3]       # rho_2(OS)
        event_avg_data[:, 6] += temp_data[:, 3]**2.
        event_avg_data[:, 7] += temp_data[:, 4]       # rho_1^2(OS)
        event_avg_data[:, 8] += temp_data[:, 4]**2.
        event_avg_data[:, 11] += temp_data[:, 6]       # rho_2(SS)
        event_avg_data[:, 12] += temp_data[:, 6]**2.
        event_avg_data[:, 13] += temp_data[:, 7]       # rho_1^2(SS)
        event_avg_data[:, 14] += temp_data[:, 7]**2.

    event_avg_data = event_avg_data/nev
    # rho_2(OS)_err
    event_avg_data[:, 6] = sqrt(event_avg_data[:, 6]
                                - event_avg_data[:, 5]**2.)/sqrt(nev)
    # rho_1^2(OS)_err
    event_avg_data[:, 8] = sqrt(event_avg_data[:, 8]
                                - event_avg_data[:, 7]**2.)/sqrt(nev)
    # rho_2(SS)_err
    event_avg_data[:, 12] = sqrt(event_avg_data[:, 12]
                                - event_avg_data[:, 11]**2.)/sqrt(nev)
    # rho_1^2(SS)_err
    event_avg_data[:, 14] = sqrt(event_avg_data[:, 14]
                                - event_avg_data[:, 13]**2.)/sqrt(nev)
    # C_2(OS)
    norm = sum(event_avg_data[:, 7])/sum(event_avg_data[:, 5])
    event_avg_data[:, 3] = norm*event_avg_data[:, 5]/event_avg_data[:, 7]
    event_avg_data[:, 4] = norm*sqrt(
          (event_avg_data[:, 6]/event_avg_data[:, 7])**2.
        + (event_avg_data[:, 5]*event_avg_data[:, 8]/event_avg_data[:, 7]**2.)**2.)
    # C_2(SS)
    norm = sum(event_avg_data[:, 13])/sum(event_avg_data[:, 11])
    event_avg_data[:, 9] = norm*event_avg_data[:, 11]/event_avg_data[:, 13]
    event_avg_data[:, 10] = norm*sqrt(
          (event_avg_data[:, 12]/event_avg_data[:, 13])**2.
        + (event_avg_data[:, 11]*event_avg_data[:, 14]/event_avg_data[:, 13]**2.)**2.)
    # Delta C_2
    event_avg_data[:, 1] = event_avg_data[:, 3] - event_avg_data[:, 9]
    event_avg_data[:, 2] = sqrt(event_avg_data[:, 4]**2.
                                + event_avg_data[:, 10]**2.)

    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ',
            header=("DeltaPhi  Delta C_2  Delta C_2_err  "
                    + "C_2(OS)  C_2(OS)_err  rho_2(OS)  rho_2(OS)_err  "
                    + "rho_1^2(OS)  rho_1^2(OS)_err  "
                    + "C_2(SS)  C_2(SS)_err  rho_2(SS)  rho_2(SS)_err  "
                    + "rho_1^2(SS)  rho_1^2(SS)_err"))
    shutil.move(file_name, avg_folder)

    file_name     = 'Correlation_function_%s_2D.dat' % ipart
    event_avg_data = loadtxt(
            path.join(file_folder_list[0],
                      'Correlation_function_9998_-9998_2D.dat'))*0.0
    event_avg_data = zeros([len(event_avg_data[:, 0]), 16])
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print("processing %s/%s ..." % (results_folder, file_name))
        temp_data = loadtxt(path.join(results_folder, file_name))
        event_avg_data[:, 0] += temp_data[:, 0]
        event_avg_data[:, 1] += temp_data[:, 1]
        event_avg_data[:, 6] += temp_data[:, 3]       # rho_2(OS)
        event_avg_data[:, 7] += temp_data[:, 3]**2.
        event_avg_data[:, 8] += temp_data[:, 4]       # rho_1^2(OS)
        event_avg_data[:, 9] += temp_data[:, 4]**2.
        event_avg_data[:, 12] += temp_data[:, 6]       # rho_2(SS)
        event_avg_data[:, 13] += temp_data[:, 6]**2.
        event_avg_data[:, 14] += temp_data[:, 7]       # rho_1^2(SS)
        event_avg_data[:, 15] += temp_data[:, 7]**2.

    event_avg_data = event_avg_data/nev
    # rho_2(OS)_err
    event_avg_data[:, 7] = sqrt(event_avg_data[:, 7]
                                - event_avg_data[:, 6]**2.)/sqrt(nev)
    # rho_1^2(OS)_err
    event_avg_data[:, 9] = sqrt(event_avg_data[:, 9]
                                - event_avg_data[:, 8]**2.)/sqrt(nev)
    # rho_2(SS)_err
    event_avg_data[:, 13] = sqrt(event_avg_data[:, 13]
                                 - event_avg_data[:, 12]**2.)/sqrt(nev)
    # rho_1^2(SS)_err
    event_avg_data[:, 15] = sqrt(event_avg_data[:, 15]
                                 - event_avg_data[:, 14]**2.)/sqrt(nev)
    # C_2(OS)
    norm = sum(event_avg_data[:, 8])/sum(event_avg_data[:, 6])
    event_avg_data[:, 4] = norm*event_avg_data[:, 6]/event_avg_data[:, 8]
    event_avg_data[:, 5] = norm*sqrt(
          (event_avg_data[:, 7]/event_avg_data[:, 8])**2.
        + (event_avg_data[:, 6]*event_avg_data[:, 9]/event_avg_data[:, 8]**2.)**2.)
    # C_2(SS)
    norm = sum(event_avg_data[:, 14])/sum(event_avg_data[:, 12])
    event_avg_data[:, 10] = norm*event_avg_data[:, 12]/event_avg_data[:, 14]
    event_avg_data[:, 11] = norm*sqrt(
          (event_avg_data[:, 13]/event_avg_data[:, 14])**2.
        + (event_avg_data[:, 12]*event_avg_data[:, 15]/event_avg_data[:, 14]**2.)**2.)
    # Delta C_2
    event_avg_data[:, 2] = event_avg_data[:, 4] - event_avg_data[:, 10]
    event_avg_data[:, 3] = sqrt(event_avg_data[:, 5]**2.
                                + event_avg_data[:, 11]**2.)

    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ',
            header=("DeltaY  DeltaPhi  Delta C_2  Delta C_2_err  "
                    + "C_2(OS)  C_2(OS)_err  rho_2(OS)  rho_2(OS)_err  "
                    + "rho_1^2(OS)  rho_1^2(OS)_err  "
                    + "C_2(SS)  C_2(SS)_err  rho_2(SS)  rho_2(SS)_err  "
                    + "rho_1^2(SS)  rho_1^2(SS)_err"))
    shutil.move(file_name, avg_folder)

print("Analysis is done.")

