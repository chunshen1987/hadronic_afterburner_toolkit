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
    num_err = zeros(len(event_avg_data[:, 0]))
    denorm_err = zeros(len(event_avg_data[:, 0]))
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        print results_folder, file_name
        temp_data = loadtxt(path.join(results_folder, file_name))
        event_avg_data += temp_data
        num_err += temp_data[:, 3]**2.
        denorm_err += temp_data[:, 4]**2.
    event_avg_data = event_avg_data/nev
    num_err = (num_err/nev - event_avg_data[:, 3]**2)/sqrt(nev)
    denorm_err = (denorm_err/nev - event_avg_data[:, 4]**2)/sqrt(nev)
    event_avg_data[:, 5] = event_avg_data[:, 3]/event_avg_data[:, 4]
    event_avg_data[:, 6] = sqrt(
        (num_err/event_avg_data[:, 4])**2
        + (event_avg_data[:, 3]*denorm_err/(event_avg_data[:, 4]**2))**2)
    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ')


print "Analysis is done."
