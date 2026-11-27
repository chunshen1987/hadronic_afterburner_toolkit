#! /usr/bin/env python
"""
     This script performs event averaging for the HBT correlation function
     calculated from event-by-event simulations

     Format from HBT correlation function file:
     col 0-2: q_out(GeV), q_side(GeV), q_long (GeV)
     col 3: numerator number of pairs in the q bin
     col 4: sum of cos(q*r) over the same-event pairs in the q bin
     col 5: denominator number of pairs in the q bin
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
import numpy as np

from centrality import multiplicity_bounds
import h5py
import shutil

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

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
    avg_folder_header = path.join(path.abspath(argv[2]),
                                  "{}_{}".format(resultsFolderName,
                                                 RapTrigLabel))
    avgResH5Filename = path.join("{}_{}_HBT.h5".format(resultsFolderName,
                                                       RapTrigLabel))
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
    avgResH5File = h5py.File("{}/{}".format(avg_folder_header,
                                            avgResH5Filename), "w")
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
    if eventStatus and file_name in event_group:
        temp_data   = event_group.get(file_name)
        temp_data   = np.asarray(temp_data)
        multiplicity = temp_data[0, 1]
        if np.isfinite(multiplicity):
            dNdyDict[event_name] = multiplicity
dNdyList = sorted(dNdyDict.values(), reverse=True)
print("Number of good events: {}".format(len(dNdyList)))



HBT_PREFIX = "HBT_correlation_function_KT_"
RAW_PREFIX = "HBT_raw_correlation_function_KT_"


def hbt_input_names(event_group):
    names = {}
    for key in event_group:
        if key.startswith(HBT_PREFIX) and key.endswith(".dat"):
            names[key] = key
    for key in event_group:
        if key.startswith(RAW_PREFIX) and key.endswith(".dat"):
            names[HBT_PREFIX + key[len(RAW_PREFIX):]] = key
    return names


def average_correlation(events, filename):
    stats = None
    raw_grid = None
    for group, names in events:
        if filename not in names:
            raise ValueError(f"{group.name} is missing {filename}")
        key = names[filename]
        table = np.asarray(group[key], dtype=float)
        if table.ndim != 2 or table.shape[1] not in (6, 8):
            raise ValueError(
                f"{group.name}/{key} has no raw pair counts. "
                "Regenerate HBT output with hbt_raw_output=1."
            )
        if not np.all(np.isfinite(table[:, 3:6])) or np.any(table[:, (3, 5)] < 0):
            raise ValueError(f"Invalid pair data in {group.name}/{key}")
        if stats is None:
            mean = np.zeros((len(table), 8))
            stats = np.zeros((6, len(table)))
        if len(table) != len(mean):
            raise ValueError(f"Inconsistent HBT grid shape in {group.name}/{key}")
        if key.startswith(RAW_PREFIX):
            if raw_grid is None:
                raw_grid = table[:, :3]
            elif not np.array_equal(raw_grid, table[:, :3]):
                raise ValueError(f"Inconsistent HBT grid in {group.name}/{key}")

        mean[:, :6] += table[:, :6]
        same, numerator, mixed = table[:, 3:6].T
        total_same, total_mixed = same.sum(), mixed.sum()
        stats[0] += numerator**2
        stats[1] += total_same**2
        stats[2] += numerator * total_same
        stats[3] += mixed**2
        stats[4] += total_mixed**2
        stats[5] += mixed * total_mixed

    nev = len(events)
    mean /= nev
    stats /= nev
    if raw_grid is not None:
        mean[:, :3] = raw_grid
    mean[:, 6:] = np.nan
    same, numerator, mixed = mean[:, 3:6].T
    total_same, total_mixed = same.sum(), mixed.sum()
    if total_same > 0 and total_mixed > 0:
        valid = (same > 0) & (mixed > 0)
        num_fraction = numerator[valid] / total_same
        den_fraction = mixed[valid] / total_mixed
        correlation = num_fraction / den_fraction
        mean[valid, 6] = correlation
        if nev > 1:
            # Propagate each normalized sum with its own covariance.
            # The last step retains the independent num/den error convention.
            var_num = (stats[0, valid] - 2 * num_fraction * stats[2, valid]
                       + num_fraction**2 * stats[1, valid]) / total_same**2
            var_den = (stats[3, valid] - 2 * den_fraction * stats[5, valid]
                       + den_fraction**2 * stats[4, valid]) / total_mixed**2
            mean[valid, 7] = np.sqrt(
                (np.maximum(var_num, 0) + correlation**2 * np.maximum(var_den, 0))
                / nev
            ) / den_fraction
    return mean


for icen in range(len(centralityCutList) - 1):
    if centralityCutList[icen+1] < centralityCutList[icen]: continue
    groupName = "{0:02.0f}-{1:02.0f}".format(centralityCutList[icen],
                                             centralityCutList[icen+1])
    cenGroup = avgResH5File.create_group(groupName)

    selected_events_list = []
    dN_dy_cut_low, dN_dy_cut_high = multiplicity_bounds(
        dNdyList, centralityCutList[icen], centralityCutList[icen+1]
    )
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
    if nev == 0:
        print("Skip ...")
        continue

    events = [(hf[name], hbt_input_names(hf[name])) for name in selected_events_list]
    filenames = sorted({filename for _, names in events for filename in names})
    if not filenames:
        raise ValueError(f"No HBT correlation datasets in centrality group {groupName}")
    for filename in filenames:
        cenGroup.create_dataset(
            filename, data=average_correlation(events, filename),
            compression="gzip", compression_opts=9,
        )

hf.close()
avgResH5File.close()
print("Analysis is done.")
