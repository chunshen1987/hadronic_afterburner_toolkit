#! /usr/bin/env python3
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
from os import path, mkdir
from glob import glob
from numpy import *
import h5py
import shutil

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

try:
    data_path = path.abspath(argv[1])
    data_name = data_path.split("/")[-1]
    results_folder_name = data_name.split(".h5")[0]
    avg_folder = path.join(path.abspath(argv[2]),
                           results_folder_name)
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
    print("Usage: {0} working_folder results_folder".format(argv[0]))
    exit(1)

hf = h5py.File(data_path, "r")
event_list = list(hf.keys())
nev = len(event_list)
dNdx_list = ['dNdtau', 'dNdx', 'dNdetas']
header_list = ["tau(fm)  dN/dtau(fm^-1)  dN/dtau_err(fm^-1)",
               (  "x(fm)  dN/dx(fm^-1)  dN/dx_err(fm^-1)  "
                + "y(fm)  dN/dy(fm^-1)  dN/dy_err(fm^-1)  "
                + "r(fm)  dN/dr(fm^-1)  dN/dr_err(fm^-1)"),
               "eta  dN/deta  dN/deta_err"]
pid = "211"

for itype, dNdx_type in enumerate(dNdx_list):
    file_name = "check_{0}_{1}.dat".format(pid, dNdx_type)
    event_group = hf.get(event_list[0])
    event_avg_data = array(event_group.get(file_name))*0.0
    for ifolder, event_name in enumerate(event_list):
        event_group = hf.get(event_name)
        print("processing {0}/{1} ...".format(event_name, file_name))
        temp_data = event_group.get(file_name)
        event_avg_data[:, :-1] += temp_data[:, :-1]
        event_avg_data[:, -1]  += temp_data[:, -1]**2.

    event_avg_data[:, -1] = sqrt(event_avg_data[:, -1])
    event_avg_data = event_avg_data/nev

    savetxt(file_name, event_avg_data, fmt='%.10e', delimiter='  ',
            header="{0}".format(header_list[itype]))

    shutil.move(file_name, avg_folder)

print("Analysis is done.")

