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
    avg_folder = path.join(path.abspath(argv[2]), working_folder.split('/')[-1])
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
    event_avg_data = zeros([len(event_avg_data[:, 0]), 3])
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print("processing %s/%s ..." % (results_folder, file_name))
        temp_data = loadtxt(path.join(results_folder, file_name))
        event_avg_data[:, 0] += temp_data[:, 0]
        event_avg_data[:, 1] += temp_data[:, 1]
        event_avg_data[:, 2] += temp_data[:, 1]**2.

    event_avg_data = event_avg_data/nev
    event_avg_data[:, 2] = sqrt(event_avg_data[:, 2]
                                - event_avg_data[:, 1]**2.)/sqrt(nev)

    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ')
    shutil.move(file_name, avg_folder)
    
    file_name = 'Balance_function_%s_Delta_phi.dat' % ipart
    event_avg_data = loadtxt(
            path.join(file_folder_list[0],
                      'Balance_function_9998_-9998_Delta_phi.dat'))*0.0
    event_avg_data = zeros([len(event_avg_data[:, 0]), 3])
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print("processing %s/%s ..." % (results_folder, file_name))
        temp_data = loadtxt(path.join(results_folder, file_name))
        event_avg_data[:, 0] += temp_data[:, 0]
        event_avg_data[:, 1] += temp_data[:, 1]
        event_avg_data[:, 2] += temp_data[:, 1]**2.

    event_avg_data = event_avg_data/nev
    event_avg_data[:, 2] = sqrt(event_avg_data[:, 2]
                                - event_avg_data[:, 1]**2.)/sqrt(nev)

    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ')
    shutil.move(file_name, avg_folder)
    
    file_name     = 'Balance_function_%s_os_2D.dat' % ipart
    file_name_err = 'Balance_function_%s_os_2D_err.dat' % ipart
    event_avg_data = loadtxt(
            path.join(file_folder_list[0],
                      'Balance_function_9998_-9998_os_2D.dat'))*0.0
    event_avg_data_err = 0.0*event_avg_data
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print("processing %s/%s ..." % (results_folder, file_name))
        temp_data = loadtxt(path.join(results_folder, file_name))
        event_avg_data     += temp_data
        event_avg_data_err += temp_data**2.

    event_avg_data     = event_avg_data/nev
    event_avg_data_err = event_avg_data_err/nev
    event_avg_data_err = sqrt((event_avg_data_err - event_avg_data**2.)/nev)

    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ')
    shutil.move(file_name, avg_folder)
    savetxt(file_name_err, event_avg_data_err, fmt='%.10e', delimiter='  ')
    shutil.move(file_name_err, avg_folder)
    
    file_name     = 'Balance_function_%s_ss_2D.dat' % ipart
    file_name_err = 'Balance_function_%s_ss_2D_err.dat' % ipart
    event_avg_data = loadtxt(
            path.join(file_folder_list[0],
                      'Balance_function_9998_-9998_ss_2D.dat'))*0.0
    event_avg_data_err = 0.0*event_avg_data
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print("processing %s/%s ..." % (results_folder, file_name))
        temp_data = loadtxt(path.join(results_folder, file_name))
        event_avg_data     += temp_data
        event_avg_data_err += temp_data**2.

    event_avg_data     = event_avg_data/nev
    event_avg_data_err = event_avg_data_err/nev
    event_avg_data_err = sqrt((event_avg_data_err - event_avg_data**2.)/nev)

    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ')
    shutil.move(file_name, avg_folder)
    savetxt(file_name_err, event_avg_data_err, fmt='%.10e', delimiter='  ')
    shutil.move(file_name_err, avg_folder)

print("Analysis is done.")

