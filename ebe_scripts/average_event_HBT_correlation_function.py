#! /usr/bin/env python
"""
     This script performs event averaging for the HBT correlation function
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
    print("Usage: average_event_HBT_correlation_function.py working_folder")
    exit(1)

file_folder_list = glob(path.join(working_folder, '*'))
results_folder_name = 'HBT_results'
KT_values = ['0_0.2', '0.2_0.4', '0.4_0.6', '0.6_0.8', '0.8_1']

nev = len(file_folder_list)
for iKT in range(len(KT_values)):
    file_name = 'HBT_correlation_function_KT_%s.dat' % KT_values[iKT]
    event_avg_data = loadtxt(path.join(file_folder_list[0],
                                       'HBT_correlation_function_KT_0_0.2.dat'))*0.0
    num = zeros(len(event_avg_data[:, 0]))
    num_err = zeros(len(event_avg_data[:, 0]))
    Npair_num = 0
    Npair_num_err = 0
    denorm = zeros(len(event_avg_data[:, 0]))
    denorm_err = zeros(len(event_avg_data[:, 0]))
    Npair_denorm = 0
    Npair_denorm_err = 0
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print "processing %s ..." % file_name
        temp_data = loadtxt(path.join(results_folder, file_name))
        num += temp_data[:, 4]
        num_err += temp_data[:, 4]**2.
        Npair_num += sum(temp_data[:, 3])
        Npair_num_err += sum(temp_data[:, 3])**2.
        denorm += temp_data[:, 5]
        denorm_err += temp_data[:, 5]**2.
        Npair_denorm += sum(temp_data[:, 5])
        Npair_denorm_err += sum(temp_data[:, 5])**2.
        event_avg_data += temp_data

    event_avg_data = event_avg_data/nev
    num = num/nev
    denorm = denorm/nev
    Npair_num = Npair_num/nev
    Npair_denorm = Npair_denorm/nev
    num_err = (num_err/nev - num**2)/sqrt(nev)
    denorm_err = (denorm_err/nev - denorm**2)/sqrt(nev)
    Npair_num_err = (Npair_num_err/nev - Npair_num**2.)/sqrt(nev)
    Npair_denorm_err = (Npair_denorm_err/nev - Npair_denorm**2.)/sqrt(nev)
    correlation = num/(Npair_num/Npair_denorm*denorm)
    correlation_err = sqrt((num_err/(Npair_num/Npair_denorm*denorm))**2. 
        + ((num*denorm_err)/((Npair_num/Npair_denorm)*denorm*denorm))**2.)
    event_avg_data[:, 4] = num
    event_avg_data[:, 5] = denorm
    event_avg_data[:, 6] = correlation
    event_avg_data[:, 7] = correlation_err
    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ')

print "Analysis is done."
