#!/usr/bin/env python

from sys import argv, exit
from os import path
from glob import glob
from numpy import *


try:
    working_folder = path.abspath(argv[1])
except(IndexError):
    print("Usage: %s working_folder" % argv[0])
    exit(1)

hbarC = 0.19733
bulk_flag = False
# convert hyper-surface from VISH2+1 format to MUSIC format
surface_file = loadtxt(path.join(working_folder, 'surface.dat'))
decdat_file = loadtxt(path.join(working_folder, 'decdat2.dat'))

n_cells = len(decdat_file[:, 0])
print("number of fluid cells: %d" % n_cells)


of = open(path.join(working_folder, "surface_MUSIC.dat"), "w")
for i in range(n_cells):
    tau_f = decdat_file[i, 0]
    x_f = surface_file[i, 2]
    y_f = surface_file[i, 3]
    eta_f = 0.0
    da0 = decdat_file[i, 1]
    da1 = decdat_file[i, 2]
    da2 = decdat_file[i, 3]
    da3 = 0.0
    v_x = decdat_file[i, 4]
    v_y = decdat_file[i, 5]
    u_tau = 1./sqrt(1. - v_x*v_x - v_y*v_y)
    u_x = u_tau*v_x
    u_y = u_tau*v_y
    u_eta = 0.0
    e_f = decdat_file[i, 6]/hbarC       # 1/fm^4
    T_f = decdat_file[i, 8]/hbarC       # 1/fm
    muB = decdat_file[i, 9]/hbarC       # 1/fm
    enthorpy_over_T = (e_f + decdat_file[i, 11]/hbarC)/T_f
    pi33 = decdat_file[i, 12]/hbarC     # 1/fm^4       
    pi00 = decdat_file[i, 13]/hbarC     # 1/fm^4
    pi01 = decdat_file[i, 14]/hbarC     # 1/fm^4
    pi02 = decdat_file[i, 15]/hbarC     # 1/fm^4
    pi03 = 0.0
    pi11 = decdat_file[i, 16]/hbarC     # 1/fm^4
    pi12 = decdat_file[i, 17]/hbarC     # 1/fm^4
    pi13 = 0.0
    pi22 = decdat_file[i, 18]/hbarC     # 1/fm^4
    pi23 = 0.0
    bulkPi = decdat_file[i, 19]/hbarC   # 1/fm^4

    of.write("%.8e  %.8e  %.8e  %.8e  " % (tau_f, x_f, y_f, eta_f))
    of.write("%.8e  %.8e  %.8e  %.8e  " % (da0, da1, da2, da3))
    of.write("%.8e  %.8e  %.8e  %.8e  " % (u_tau, u_x, u_y, u_eta))
    of.write("%.8e  %.8e  %.8e  %.8e  " % (e_f, T_f, muB, enthorpy_over_T))
    of.write("%.8e  %.8e  %.8e  %.8e  " % (pi00, pi01, pi02, pi03))
    of.write("%.8e  %.8e  %.8e  " % (pi11, pi12, pi13))
    of.write("%.8e  %.8e  %.8e" % (pi22, pi23, pi33))
    if (bulk_flag):
        of.write("  %.8e" % bulkPi)
    of.write("\n")
of.close()

