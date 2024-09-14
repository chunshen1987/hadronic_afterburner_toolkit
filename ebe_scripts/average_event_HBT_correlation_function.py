#! /usr/bin/env python
"""
     This script performs event averaging for the HBT correlation function
     calculated from event-by-event simulations

     Format from HBT correlation function file:
     (col 0-2: q_out(GeV), q_side(GeV), q_long (GeV) ecoOutput == 0)
     col 0: numerator <cos(q*r)> in the q bin
     col 1: denominator number of pairs in the q bin
"""

from sys import argv, exit
from os import path, makedirs
from glob import glob
import numpy as np

try:
    working_folder = path.abspath(argv[1])
    avg_folder = path.abspath(argv[2])
    makedirs(avg_folder, exist_ok=True)
except IndexError:
    print(f"Usage: {argv[0]} working_folder results_folder")
    exit(1)


file_folder_list = glob(path.join(working_folder, 'UrQMD*'))
results_folder_name = 'UrQMD_results'
hbtFileList = [
    ifile.split("/")[-1] for ifile in glob(
                path.join(file_folder_list[0], results_folder_name, "HBT*"))]

nev = len(file_folder_list)
for hbtFile_i in hbtFileList:
    event_avg_data = 0.*np.loadtxt(
        path.join(file_folder_list[0], results_folder_name, hbtFile_i))
    for folder_i in file_folder_list:
        results_folder = path.abspath(path.join(folder_i, results_folder_name))
        print(f"processing {results_folder}/{hbtFile_i} ...")
        temp_data = np.loadtxt(path.join(results_folder, hbtFile_i))
        event_avg_data += temp_data

    event_avg_data = event_avg_data/nev

    np.savetxt(path.join(avg_folder, hbtFile_i), event_avg_data,
               fmt='%.10e', delimiter='  ')

print("Analysis is done.")

