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

KT_values = ['0.15_0.25', '0.25_0.35', '0.35_0.45', '0.45_0.55']

Reg_centrality_cut_list = [0., 5., 10., 20., 30., 40., 50.,
                           60., 70., 80., 90., 100.]
PHOBOS_cen_list = [0., 6., 15., 25., 35., 45., 55.]  # PHOBOS AuAu 200
SPS_cen_list    = [5., 12.5, 23.5, 33.5, 43.5]       # SPS PbPb
PHENIX_cen_list = [0., 20., 40., 60., 88.]           # PHENIX dAu
STAR_cen_list   = [0., 10., 40., 80]                 # STAR v1
ALICE_pp_list   = [0., 100., 0., 1., 5.,
                   0., 5., 10., 15, 20., 30., 40., 50., 70., 100.]
centralityCutList = Reg_centrality_cut_list
dNcutList = []    # pre-defined Nch cut if simulation is not minimum bias

RapidityTrigger = 0  # 0: mid-rapidity [-0.5, 0.5]
                     # 1: PHENIX BBC trigger [-3.9, -3.1]
                     # 2: ALICE V0A trigger [-5.1, -2.8]
                     # 3: ATLAS forward trigger [-4.9, -3.1]

RapTrigLabel = "CL1"
if RapidityTrigger == 1:
    RapTrigLabel = "BBC"
elif RapidityTrigger == 2:
    RapTrigLabel = "V0A"
elif RapidityTrigger == 3:
    RapTrigLabel = "ATLASForward"

try:
    data_path = path.abspath(argv[1])
    data_name = data_path.split("/")[-1]
    data_path = "/".join(data_path.split("/")[0:-1])
    resultsFolderName = data_name.split(".h5")[0]
    avg_folder_header = path.join(
        path.abspath(argv[2]), "{}_{}".format(resultsFolderName, RapTrigLabel))
    print("output folder: %s" % avg_folder_header)
    if path.isdir(avg_folder_header):
        print("folder %s already exists!" % avg_folder_header)
        var = input("do you want to delete it? [y/N]")
        if 'y' in var.lower():
            shutil.rmtree(avg_folder_header)
            mkdir(avg_folder_header)
        else:
            print("Continue analysis in {} ...".format(avg_folder_header))
    else:
        mkdir(avg_folder_header)
    paraFiles = glob(path.join(data_path, "parameter*"))
    for iFile in paraFiles:
        shutil.copy(iFile, avg_folder_header)
except IndexError:
    print("Usage: {} working_folder results_folder".format(argv[0]))
    exit(1)


def check_an_event_is_good(h5_event):
    """This function checks the given event contains all required files"""
    required_files_list = [
        'particle_9999_vndata_eta_-0.5_0.5.dat',
        'particle_211_vndata_diff_y_-0.5_0.5.dat',
        'particle_321_vndata_diff_y_-0.5_0.5.dat',
        'particle_2212_vndata_diff_y_-0.5_0.5.dat',
        'particle_3122_vndata_diff_y_-0.5_0.5.dat',
        'particle_-3122_vndata_diff_y_-0.5_0.5.dat',
    ]
    event_file_list = list(h5_event.keys())
    for ifile in required_files_list:
        if ifile not in event_file_list:
            print("event {} is bad, missing {} ...".format(h5_event.name,
                                                           ifile))
            return False
    return True

hf = h5py.File(argv[1], "r")
event_list = list(hf.keys())
print("total number of events: {}".format(len(event_list)))
nev = len(event_list)

dNdyDict = {}
dNdyList = []
for ifolder, event_name in enumerate(event_list):
    file_name = "particle_9999_vndata_eta_-0.5_0.5.dat"
    if RapidityTrigger == 1:      # PHENIX BBC Trigger
        file_name = "particle_9999_vndata_eta_-3.9_-3.1.dat"
    elif RapidityTrigger == 2:    # ALICE V0A Trigger
        file_name = "particle_9999_vndata_eta_-5.1_-2.8.dat"
    elif RapidityTrigger == 3:    # ATLAS forward Trigger
        file_name = "particle_9999_vndata_eta_-4.9_-3.1.dat"
    event_group = hf.get(event_name)
    eventStatus = check_an_event_is_good(event_group)
    if eventStatus:
        temp_data   = event_group.get(file_name)
        temp_data   = nan_to_num(temp_data)
        dNdyDict[event_name] = temp_data[0, 1]
dNdyList = -sort(-array(list(dNdyDict.values())))
print("Number of good events: {}".format(len(dNdyList)))


for icen in range(len(centralityCutList) - 1):
    if centralityCutList[icen+1] < centralityCutList[icen]: continue
    avg_folder = path.join(
        avg_folder_header, "{0:02.0f}-{1:02.0f}".format(
            centralityCutList[icen], centralityCutList[icen+1])
    )

    if path.isdir(avg_folder):
        print("{} already exists, skipped ...".format(avg_folder))
        continue
    else:
        mkdir(avg_folder)

    selected_events_list = []
    dN_dy_cut_high = dNdyList[
        int(len(dNdyList)*centralityCutList[icen]/100.)
    ]
    dN_dy_cut_low  = dNdyList[
        min(len(dNdyList)-1,
            int(len(dNdyList)*centralityCutList[icen+1]/100.))
    ]
    if len(dNcutList) == len(centralityCutList):
        dN_dy_cut_high = dNcutList[icen]
        dN_dy_cut_low = dNcutList[icen+1]

    for event_name in dNdyDict.keys():
        if (dNdyDict[event_name] > dN_dy_cut_low
            and dNdyDict[event_name] <= dN_dy_cut_high):
            selected_events_list.append(event_name)

    nev = len(selected_events_list)
    print("analysis {}%-{}% nev = {}...".format(
            centralityCutList[icen], centralityCutList[icen+1], nev))
    print("dNdy: {0:.2f} - {1:.2f}".format(dN_dy_cut_low, dN_dy_cut_high))
    if nev == 0:
        print("Skip ...")
        continue

    for iKT in range(len(KT_values)):
        file_name = 'HBT_correlation_function_KT_%s.dat' % KT_values[iKT]
        event_group = hf.get(event_list[0])
        event_avg_data = array(event_group.get(file_name))*0.0
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

        for ifolder, event_name in enumerate(selected_events_list):
            event_group = hf.get(event_name)
            #print("processing {0}/{1} ...".format(event_name, file_name))
            temp_data = event_group.get(file_name)
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

        correlation = num/(Npair_num/Npair_denorm*denorm)

        err_numerator = ((num/Npair_num*sqrt(sigma_num/(num**2.)
                          + sigma_Npair_num/(Npair_num**2.)
                         - 2.*sigma_num_Npair_num/(num*Npair_num))/sqrt(nev)))
        err_denormator = ((denorm/Npair_denorm*sqrt(sigma_denorm/(denorm**2.)
                           + sigma_Npair_denorm/(Npair_denorm**2.)
                          - 2.*sigma_denorm_Npair_denorm/(denorm*Npair_denorm))
                          /sqrt(nev)))

        correlation_err = sqrt((err_numerator/(denorm/Npair_denorm))**2.
                               + ((num/Npair_num)*err_denormator
                                   /((denorm/Npair_denorm)**2.))**2.)
        event_avg_data[:, 4] = num
        event_avg_data[:, 5] = denorm
        event_avg_data[:, 6] = correlation
        event_avg_data[:, 7] = correlation_err
        savetxt(path.join(avg_folder, file_name), event_avg_data,
                fmt='%.10e', delimiter='  ')

print("Analysis is done.")

