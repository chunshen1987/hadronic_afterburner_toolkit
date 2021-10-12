#! /usr/bin/env python
"""
     This script performs event averaging for the HBT correlation function
     calculated from event-by-event simulations

     Format from HBT correlation function file:
     col 0-2: q_out(GeV), q_side(GeV), q_long (GeV)
     col 3: numerator number of pairs in the q bin
     col 4: numerator <cos(q*r)> in the q bin
     col 5: denominator number of pairs in the q bin
"""

from sys import argv, exit
from os import path, makedirs
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
    avg_folder = path.abspath(argv[2])
    makedirs(avg_folder, exist_ok=True)
except IndexError:
    print("Usage: average_event_HBT_correlation_function.py working_folder results_folder")
    exit(1)

file_folder_list = glob(path.join(working_folder, 'UrQMD*'))
results_folder_name = 'UrQMD_results'
KT_values = ['0_0.2', '0.2_0.4', '0.4_0.6', '0.6_0.8', '0.8_1']

nev = len(file_folder_list)
for iKT in range(len(KT_values)):
    file_name = 'HBT_correlation_function_KT_%s.dat' % KT_values[iKT]
    event_avg_data = loadtxt(path.join(
        file_folder_list[0], results_folder_name,
        'HBT_correlation_function_KT_0_0.2.dat'))*0.0
    num = zeros(len(event_avg_data[:, 0]))
    sigma_num = zeros(len(event_avg_data[:, 0]))
    Npair_num = 0
    sigma_Npair_num = 0
    sigma_num_Npair_num = zeros(len(event_avg_data[:, 0]))
    denorm = zeros(len(event_avg_data[:, 0]))
    sigma_denorm = zeros(len(event_avg_data[:, 0]))
    sigma_denorm_Npair_denorm = zeros(len(event_avg_data[:, 0]))
    Npair_denorm = 0
    sigma_Npair_denorm = 0

    for ifolder in range(nev):
        results_folder = path.abspath(
                path.join(file_folder_list[ifolder], results_folder_name))
        print("processing %s/%s ..." % (results_folder, file_name))
        temp_data = loadtxt(path.join(results_folder, file_name))
        num += temp_data[:, 4]
        sigma_num += temp_data[:, 4]**2.
        Npair_num += sum(temp_data[:, 3])
        sigma_Npair_num += sum(temp_data[:, 3])**2.
        sigma_num_Npair_num += temp_data[:, 4]*sum(temp_data[:, 3])
        denorm += temp_data[:, 5]
        sigma_denorm += temp_data[:, 5]**2.
        Npair_denorm += sum(temp_data[:, 5])
        sigma_Npair_denorm += sum(temp_data[:, 5])**2.
        sigma_denorm_Npair_denorm += temp_data[:, 5]*sum(temp_data[:, 5])
        event_avg_data += temp_data

    event_avg_data = event_avg_data/nev
    num = num/nev
    denorm = denorm/nev
    Npair_num = Npair_num/nev
    Npair_denorm = Npair_denorm/nev
    sigma_num = sigma_num/nev - num**2
    sigma_denorm = sigma_denorm/nev - denorm**2
    sigma_Npair_num = sigma_Npair_num/nev - Npair_num**2.
    sigma_Npair_denorm = sigma_Npair_denorm/nev - Npair_denorm**2.
    sigma_num_Npair_num = sigma_num_Npair_num/nev - num*Npair_num
    sigma_denorm_Npair_denorm = (sigma_denorm_Npair_denorm/nev
                                 - denorm*Npair_denorm)

    correlation = nan_to_num(num/(Npair_num/Npair_denorm*denorm))

    err_numerator = nan_to_num(
        num/Npair_num*sqrt(sigma_num/(num**2.) + sigma_Npair_num/(Npair_num**2.)
        - 2.*sigma_num_Npair_num/(num*Npair_num))/sqrt(nev))
    err_denormator = nan_to_num(
        (denorm/Npair_denorm*sqrt(sigma_denorm/(denorm**2.)
            + sigma_Npair_denorm/(Npair_denorm**2.)
        - 2.*sigma_denorm_Npair_denorm/(denorm*Npair_denorm))/sqrt(nev)))

    correlation_err = nan_to_num(
        sqrt((err_numerator/(denorm/Npair_denorm))**2.
        + ((num/Npair_num)*err_denormator/((denorm/Npair_denorm)**2.))**2.))
    event_avg_data[:, 4] = num
    event_avg_data[:, 5] = denorm
    event_avg_data[:, 6] = correlation
    event_avg_data[:, 7] = correlation_err
    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ')

    shutil.move(file_name, avg_folder)

print("Analysis is done.")

