#!/usr/bin/env python

from os import path
from numpy import loadtxt

collision_system = "PbPb2760"
centrality_list = ['0-5', '5-10', '10-20', '20-30', '30-40', '40-50']
centrality_mid = [2.5, 7.5, 15, 25, 35, 45]
folder_basename = "%s_IPG_C%s_Tsw_145"
file_name = "charged_hadron_event_plane_correlation_ATLAS.dat"

output_file = open("%s_charged_hadron_event_plane_correlation_ATLAS.dat"
                   % collision_system, "w")
# write the header
output_file.write("# centrality 4(24) 4(24)_err 6(23) 6(23)_err "
                  "6(26) 6(26)_err 6(36) 6(36)err (235) (235)_err "
                  "(246) (246)_err (234) (234)_err\n")
for i, centrality_text in enumerate(centrality_list):
    folder_name = folder_basename % (collision_system, centrality_text)
    data = open(path.join('.', folder_name, file_name), "r")
    data.readline()
    output_file.write("%g  " % (centrality_mid[i]))
    for line in data.readlines():
        value = float(line.split(" ")[2])
        value_err = float(line.split(" ")[4].split("\n")[0])
        output_file.write("%.6e  %.6e  " % (value, value_err))
    data.close()
    output_file.write("\n")
output_file.close()

file_name = "charged_hadron_event_plane_correlation_ALICE.dat"
output_file = open("%s_charged_hadron_event_plane_correlation_ALICE.dat"
                   % collision_system, "w")
# write the header
output_file.write("# centrality 4(24) 4(24)_err 6(23) 6(23)_err "
                  "6(26) 6(26)_err 6(36) 6(36)err (235) (235)_err "
                  "(246) (246)_err (234) (234)_err\n")
for i, centrality_text in enumerate(centrality_list):
    folder_name = folder_basename % (collision_system, centrality_text)
    data = open(path.join('.', folder_name, file_name), "r")
    data.readline()
    output_file.write("%g  " % (centrality_mid[i]))
    for line in data.readlines():
        value = float(line.split(" ")[2])
        value_err = float(line.split(" ")[4].split("\n")[0])
        output_file.write("%.6e  %.6e  " % (value, value_err))
    data.close()
    output_file.write("\n")
output_file.close()
