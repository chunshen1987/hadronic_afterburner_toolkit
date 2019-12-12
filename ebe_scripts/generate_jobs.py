#!/usr/bin/env python

import sys
from os import path, mkdir
import shutil
from glob import glob
import subprocess
import random

def write_script_header(cluster, script, event_id, walltime, working_folder):
    if cluster == "nersc":
        script.write(
"""#!/bin/bash -l
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -J UrQMD_%s
#SBATCH -t %s
#SBATCH -L SCRATCH
#SBATCH -C haswell
""" % (event_id, walltime))
    elif cluster == "guillimin":
        script.write(
"""#!/usr/bin/env bash
#PBS -N UrQMD_%s
#PBS -l nodes=1:ppn=1
#PBS -l walltime=%s
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ad
#PBS -q sw
#PBS -d %s
""" % (event_id, walltime, working_folder))
    elif cluster == "McGill":
        script.write(
"""#!/usr/bin/env bash
#PBS -N UrQMD_%s
#PBS -l nodes=1:ppn=1:irulan
#PBS -l walltime=%s
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -d %s
""" % (event_id, walltime, working_folder))
    else:
        print("Error: unrecoginzed cluster name :", cluster)
        print("Available options: nersc, guillimin, McGill")
        exit(1)


def write_analysis_spectra_and_vn_commands(script, after_burner_type):
    pid_particle_list = ['211', '-211', '321', '-321', '2212', '-2212',
                         '3122', '-3122', '3312', '-3312', '3334', '-3334',
                         '333']
    charged_particle_list = ['9999', '9998', '-9998']
    #pid_particle_list = []
    #charged_particle_list = ['9999']

    read_in_mode = 2
    if after_burner_type == "JAM":
        read_in_mode = 5
    if after_burner_type == "OSCAR":
        read_in_mode = 0

    for ipart in charged_particle_list:
        script.write(
"""
    # charged hadrons
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=-0.1 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.1 rap_max=1.0 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=0.5 rap_max=2.0 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=-0.5 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0 compute_correlation=1 flag_charge_dependence=1 pT_min=0.2 pT_max=2.0 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=2.0 compute_correlation=1 flag_charge_dependence=1 pT_min=0.2 pT_max=2.0 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-2.0 rap_max=2.0 >> ../output.log
""".format(read_in_mode, ipart))
    for ipart in pid_particle_list:
        script.write(
"""
    #./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 >> ../output.log
    ./hadronic_afterburner_tools.e run_mode=0 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 >> ../output.log
""".format(read_in_mode, ipart))


def write_analysis_particle_distrubtion_commands(script, after_burner_type):
    pid_particle_list = ['211', '-211', '321', '-321', '2212', '-2212',
                         '3122', '-3122']
    charged_particle_list = ['9997', '-9997', '9998', '-9998']

    read_in_mode = 2
    if after_burner_type == "JAM":
        read_in_mode = 5
    if after_burner_type == "OSCAR":
        read_in_mode = 0

    for ipart in pid_particle_list:
        script.write(
"""
    ./hadronic_afterburner_tools.e run_mode=2 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=0 >> output.log
    ./hadronic_afterburner_tools.e run_mode=2 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=1 rap_type=1 >> output.log
""".format(read_in_mode, ipart))
        if "-" not in ipart:
            script.write(
"""
    ./hadronic_afterburner_tools.e run_mode=2 read_in_mode={0} particle_monval={1} distinguish_isospin=1 rap_type=0 net_particle_flag=1 >> output.log
    ./hadronic_afterburner_tools.e run_mode=2 read_in_mode={0} particle_monval={1} distinguish_isospin=1 rap_type=1 net_particle_flag=1 >> output.log
""".format(read_in_mode, ipart))
    for ipart in charged_particle_list:
        script.write(
"""
    ./hadronic_afterburner_tools.e run_mode=2 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 >> output.log
""".format(read_in_mode, ipart))
        if "-" not in ipart:
            script.write(
"""
    ./hadronic_afterburner_tools.e run_mode=2 read_in_mode={0} particle_monval={1} resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 net_particle_flag=1 >> output.log
""".format(read_in_mode, ipart))
    script.write(
"""
    ./hadronic_afterburner_tools.e run_mode=2 read_in_mode={0} particle_monval=9999 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 >> output.log
""".format(read_in_mode))


def generate_script(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '10:00:00'
    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir UrQMD_results
for iev in `ls OSCAR_events`
do
    cd osc2u
    ./osc2u.e < ../OSCAR_events/$iev
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh
    mv particle_list.dat ../UrQMD_results/particle_list_`echo $iev | cut -f 2 -d _`
    cd ..
done
""")
    script.close()


def generate_script_JAM(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '10:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/irulan/chun/JAM/JAM_lib/lib
mkdir JAM_results
for iev in `ls OSCAR_events`
do
    eventid=`echo $iev | cut -f 2 -d "_" | cut -f 1 -d "."`
    cd JAM
    mv ../OSCAR_events/$iev ./OSCAR.DAT
    rm -fr phase.dat
    ./jamgo
    mv phase.dat ../JAM_results/particle_list_$eventid.dat
    mv OSCAR.DAT ../OSCAR_events/OSCAR_$eventid.dat
    cd ..
done
""")
    script.close()


def generate_script_iSS(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '35:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir UrQMD_results
mkdir spvn_results
for iev in `ls hydro_events --color=none | grep "surface"`
do
    event_id=`echo $iev | rev |  cut -f 1 -d _ | rev | cut -f 1 -d .`
    cd iSS
    if [ -d "results" ]; then
        rm -fr results
    fi
    mkdir results
    mv ../hydro_events/$iev results/surface.dat
    cp ../hydro_events/music_input_event_$event_id results/music_input
    ./iSS.e >> ../output.log
    mv results/surface.dat ../hydro_events/$iev
    #rm -fr results/sample*
    # turn on global momentum conservation
    ./correct_momentum_conservation.py OSCAR.DAT
    mv OSCAR_w_GMC.DAT OSCAR.DAT
    cd ../osc2u
    ./osc2u.e < ../iSS/OSCAR.DAT >> ../output.log
    mv fort.14 ../urqmd/OSCAR.input
    cd ../urqmd
    ./runqmd.sh >> ../output.log
    mv particle_list.dat ../UrQMD_results/particle_list_$event_id.dat
    #mv ../iSS/OSCAR.DAT ../UrQMD_results/OSCAR_$event_id.dat
    rm -fr ../iSS/OSCAR.DAT
    rm -fr OSCAR.input
    cd ..
    ./hadronic_afterburner_toolkit/convert_to_binary.e UrQMD_results/particle_list_$event_id.dat
    rm -fr UrQMD_results/particle_list_$event_id.dat
    
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../UrQMD_results/particle_list_$event_id.gz results/particle_list.dat
""")
    write_analysis_spectra_and_vn_commands(script, "UrQMD")
    script.write(
"""
    mv results/particle_list.dat ../UrQMD_results/particle_list_$event_id.gz
    mv results ../spvn_results/event_$event_id
    cd ..
done
""")
    script.close()


def generate_script_iS(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '3:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir spvn_results
for iev in `ls hydro_events --color=none | grep "surface"`
do
    event_id=`echo $iev | rev | cut -f 1 -d _ | rev | cut -f 1 -d .`
    cd iS
    if [ -d "results" ]; then
        rm -fr results
    fi
    mkdir results
    mv ../hydro_events/$iev results/surface.dat
    cp ../hydro_events/music_input_event_$event_id results/music_input
    ./iS_withResonance.sh >> ../output.log
    mv results/surface.dat ../hydro_events/$iev
    mv results/ ../spvn_results/event_$event_id
    cd ..
done
""")
    script.close()

def generate_script_HBT(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '20:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir HBT_results
for iev in `ls UrQMD_events | grep "particle_list"`
do
    eventid=`echo $iev | rev | cut -f 1 -d _ | rev | cut -f 1 -d .`
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../UrQMD_events/$iev results/particle_list.dat
    mv ../UrQMD_events/mixed_event_$eventid.dat results/particle_list_mixed_event.dat
    ./hadronic_afterburner_tools.e read_in_mode=2 run_mode=1 resonance_feed_down_flag=0 > output.log
    mv results/particle_list.dat ../UrQMD_events/$iev
    mv results/particle_list_mixed_event.dat ../UrQMD_events/mixed_event_$eventid.dat
    mv results ../HBT_results/event_$eventid
    cd ..
done
""")
    script.close()


def generate_script_HBT_with_JAM(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '30:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir HBT_results
for iev in `ls JAM_events | grep "particle_list"`
do
    eventid=`echo $iev | rev | cut -f 1 -d _ | rev | cut -f 1 -d .`
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../JAM_events/$iev results/particle_list.dat
    mv ../JAM_events/mixed_event_$eventid.dat results/particle_list_mixed_event.dat
    ./hadronic_afterburner_tools.e run_mode=1 read_in_mode=5 resonance_feed_down_flag=0 > output.log
    mv results/particle_list.dat ../JAM_events/$iev
    mv results/particle_list_mixed_event.dat ../JAM_events/mixed_event_$eventid.dat
    mv results ../HBT_results/event_$eventid
    cd ..
done
""")
    script.close()


def generate_script_balance_function(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '01:00:00'

    particle_a_list = ['9998']
    particle_b_list = ['-9998']

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir BalanceFunction_results
for iev in `ls UrQMD_events | grep "particle_list"`
do
    eventid=`echo $iev | rev | cut -f 1 -d _ | rev | cut -f 1 -d .`
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../UrQMD_events/$iev results/particle_list.dat
    mv ../UrQMD_events/mixed_event_$eventid.dat results/particle_list_mixed_event.dat
""")
    for ipart in range(len(particle_a_list)):
        script.write(
"""
    ./hadronic_afterburner_tools.e read_in_mode=2 run_mode=3 resonance_feed_down_flag=0 distinguish_isospin=0 rap_type=0 rap_min=-1.0 rap_max=1.0 particle_alpha={0} particle_beta={1} BpT_min=0.2 BpT_max=3.0 > output.log
""".format(particle_a_list[ipart], particle_b_list[ipart]))
    script.write(
"""
    mv results/particle_list.dat ../UrQMD_events/$iev
    mv results/particle_list_mixed_event.dat ../UrQMD_events/mixed_event_$eventid.dat
    mv results ../BalanceFunction_results/event_$eventid
    cd ..
done
""")
    script.close()

def generate_script_spectra_and_vn(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '1:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir spvn_results
for iev in `ls UrQMD_events | grep "particle_list"`
do
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../UrQMD_events/$iev results/particle_list.dat
""")
    write_analysis_spectra_and_vn_commands(script, "UrQMD")
    script.write(
"""
    mv results/particle_list.dat ../UrQMD_events/$iev
    mv results ../spvn_results/event_`echo $iev | rev | cut -f 1 -d _ | rev | cut -f 1 -d .`
    cd ..
done
""")
    script.close()


def generate_script_particle_yield_distribution(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '1:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir spvn_results
for iev in `ls UrQMD_events | grep "particle_list"`
do
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../UrQMD_events/$iev results/particle_list.dat
""")
    write_analysis_particle_distrubtion_commands(script, "UrQMD")
    script.write(
"""
    mv results/particle_list.dat ../UrQMD_events/$iev
    mv results ../spvn_results/event_`echo $iev | rev | cut -f 1 -d _ | rev | cut -f 1 -d .`
    cd ..
done
""")
    script.close()


def generate_script_particle_yield_distribution_with_OSCAR(cluster_name,
                                                           folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '1:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir spvn_results
for iev in `ls OSCAR_events`
do
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../OSCAR_events/$iev results/OSCAR.DAT
""")
    write_analysis_particle_distrubtion_commands(script, "OSCAR")
    script.write(
"""
    mv results/OSCAR.DAT ../OSCAR_events/$iev
    mv results ../spvn_results/event_`echo $iev | cut -f 2 -d _ | cut -f 1 -d .`
    cd ..
done
""")
    script.close()


def generate_script_spectra_and_vn_with_JAM(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '3:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir spvn_results
for iev in `ls JAM_events | grep "particle_list"`
do
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../JAM_events/$iev results/particle_list.dat
""")
    write_analysis_spectra_and_vn_commands(script, "JAM")
    script.write(
"""
    mv results/particle_list.dat ../JAM_events/$iev
    mv results ../spvn_results/event_`echo $iev | rev | cut -f 1 -d _ | rev | cut -f 1 -d .`
    cd ..
done
""")
    script.close()


def generate_script_HBT_with_OSCAR(cluster_name, folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    event_id = working_folder.split('/')[-1]
    walltime = '35:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    write_script_header(cluster_name, script, event_id, walltime,
                        working_folder)
    script.write(
"""
mkdir HBT_results
for iev in `ls OSCAR_events | grep "OSCAR"`
do
    eventid=`echo $iev | cut -f 2 -d _ | cut -f 1 -d .`
    cd hadronic_afterburner_toolkit
    rm -fr results
    mkdir results
    mv ../OSCAR_events/$iev results/OSCAR.DAT
    mv ../OSCAR_events/mixed_event_$eventid.dat results/OSCAR_mixed_event.DAT
    ./hadronic_afterburner_tools.e read_in_mode=0 run_mode=1 resonance_feed_down_flag=1 > output.log
    mv results/OSCAR.DAT ../OSCAR_events/$iev
    mv results/OSCAR_mixed_event.DAT ../OSCAR_events/mixed_event_$eventid.dat
    mv results ../HBT_results/event_$eventid
    cd ..
done
""")
    script.close()


def copy_UrQMD_events(number_of_cores, input_folder, working_folder):
    events_list = glob('%s/particle_list_*.dat' % input_folder)
    if events_list == []:
        events_list = glob('%s/particle_list_*.gz' % input_folder)
        if events_list == []:
            print("Error: can not find UrQMD events, events_list is empty! ",
                  events_list)
        else:
            print("Linking zipped binary UrQMD events, ",
                  "make sure read_in_mode is set to 2~")

    for iev in range(len(events_list)):
        folder_id = iev % number_of_cores
        filename = events_list[iev].split('/')[-1].split('.')[0]
        event_id = filename.split('_')[-1]
        folder_path = path.join(working_folder, 'event_%d' % folder_id, 
                                'UrQMD_events', '%s.dat' % filename)
        bashCommand = "ln -s %s %s" % (
            path.abspath(events_list[iev]), folder_path)
        subprocess.Popen(bashCommand, stdout = subprocess.PIPE, shell=True)
        mixed_id = random.randint(0, len(events_list)-1)
        filename_mixed = events_list[mixed_id].split('/')[-1].split('.')[0]
        mixed_event_id = filename_mixed.split('_')[-1]
        while (mixed_event_id == iev):
            mixed_id = random.randint(0, len(events_list)-1)
            filename_mixed = events_list[mixed_id].split('/')[-1].split('.')[0]
            mixed_event_id = filename_mixed.split('_')[-1]
        folder_path = path.join(
            working_folder, 'event_%d' % folder_id, 
            'UrQMD_events', 'mixed_event_%s.dat' % event_id)
        bashCommand = "ln -s %s %s" % (
            path.abspath(events_list[mixed_id]), folder_path)
        subprocess.Popen(bashCommand, stdout = subprocess.PIPE, shell=True)


def copy_JAM_events(number_of_cores, input_folder, working_folder):
    events_list = glob('%s/particle_list_*.dat' % input_folder)
    for iev in range(len(events_list)):
        folder_id = iev % number_of_cores
        filename = events_list[iev].split('/')[-1].split('.')[0]
        event_id = filename.split('_')[-1]
        folder_path = path.join(working_folder, 'event_%d' % folder_id, 
                                'JAM_events', '%s.dat' % filename)
        bashCommand = "ln -s %s %s" % (
            path.abspath(events_list[iev]), folder_path)
        subprocess.Popen(bashCommand, stdout = subprocess.PIPE, shell=True)
        mixed_id = random.randint(0, len(events_list)-1)
        filename_mixed = events_list[mixed_id].split('/')[-1].split('.')[0]
        mixed_event_id = filename_mixed.split('_')[-1]
        while (mixed_event_id == iev):
            mixed_id = random.randint(0, len(events_list)-1)
            filename_mixed = events_list[mixed_id].split('/')[-1].split('.')[0]
            mixed_event_id = filename_mixed.split('_')[-1]
        folder_path = path.join(working_folder, 'event_%d' % folder_id, 
                                'JAM_events', 'mixed_event_%s.dat' % event_id)
        bashCommand = "ln -s %s %s" % (path.abspath(events_list[mixed_id]),
                                       folder_path)
        subprocess.Popen(bashCommand, stdout = subprocess.PIPE, shell=True)


def generate_event_folder_UrQMD(cluster_name, working_folder, event_id, mode):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)

    if mode == 2:
        # calculate HBT correlation with OSCAR outputs
        mkdir(path.join(event_folder, 'OSCAR_events'))
        generate_script_HBT_with_OSCAR(cluster_name, event_folder)
    elif mode == 3:
        # calculate HBT correlation with UrQMD outputs
        mkdir(path.join(event_folder, 'UrQMD_events'))
        generate_script_HBT(cluster_name, event_folder)
    elif mode == 4:
        # calculate HBT correlation with UrQMD outputs
        mkdir(path.join(event_folder, 'UrQMD_events'))
        generate_script_spectra_and_vn(cluster_name, event_folder)
    elif mode == 8:
        # collect event-by-event particle distribution
        mkdir(path.join(event_folder, 'UrQMD_events'))
        generate_script_particle_yield_distribution(cluster_name, event_folder)
    elif mode == 9:
        # calculate event-by-event particle distribution with OSCAR outputs
        mkdir(path.join(event_folder, 'OSCAR_events'))
        generate_script_particle_yield_distribution_with_OSCAR(cluster_name,
                                                               event_folder)
    elif mode == 10:
        # calculate balance function correlation with UrQMD outputs
        mkdir(path.join(event_folder, 'UrQMD_events'))
        generate_script_balance_function(cluster_name, event_folder)

    shutil.copytree('codes/hadronic_afterburner_toolkit', 
                    path.join(path.abspath(event_folder), 
                    'hadronic_afterburner_toolkit'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath(path.join('codes', 'hadronic_afterburner_toolkit_code',
                               'hadronic_afterburner_tools.e')),
        path.join(path.abspath(event_folder), "hadronic_afterburner_toolkit",
                  "hadronic_afterburner_tools.e")), shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/hadronic_afterburner_toolkit_code/EOS'),
        path.join(path.abspath(event_folder),
                  "hadronic_afterburner_toolkit/EOS")), shell=True)


def generate_event_folder_JAM(cluster_name, working_folder, event_id, mode):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)

    if mode == 5:
        # run JAM with OSCAR files
        mkdir(path.join(event_folder, 'OSCAR_events'))
        generate_script_JAM(cluster_name, event_folder)
        shutil.copytree('codes/JAM', 
                        path.join(path.abspath(event_folder), 'JAM'))
    elif mode == 6:
        # collect particle spectra and vn with JAM outputs
        mkdir(path.join(event_folder, 'JAM_events'))
        generate_script_spectra_and_vn_with_JAM(cluster_name, event_folder)
        shutil.copytree('codes/hadronic_afterburner_toolkit', 
                        path.join(path.abspath(event_folder), 
                        'hadronic_afterburner_toolkit'))
    elif mode == 7:
        # calculate HBT correlation with JAM outputs
        mkdir(path.join(event_folder, 'JAM_events'))
        generate_script_HBT_with_JAM(cluster_name, event_folder)
        shutil.copytree('codes/hadronic_afterburner_toolkit', 
                        path.join(path.abspath(event_folder), 
                        'hadronic_afterburner_toolkit'))


def generate_event_folder(cluster_name, working_folder, event_id):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    mkdir(path.join(event_folder, 'OSCAR_events'))
    generate_script(cluster_name, event_folder)
    shutil.copytree('codes/osc2u', 
                    path.join(path.abspath(event_folder), 'osc2u'))
    shutil.copytree('codes/urqmd', 
                    path.join(path.abspath(event_folder), 'urqmd'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/urqmd_code/urqmd/urqmd.e'),
        path.join(path.abspath(event_folder), "urqmd/urqmd.e")), shell=True)

def copy_OSCAR_events(number_of_cores, input_folder, working_folder):
    events_list = glob('%s/*.dat' % input_folder)
    for iev in range(len(events_list)):
        folder_id = iev % number_of_cores
        filename = events_list[iev].split('/')[-1].split('.')[0] 
        event_id = filename.split('_')[-1]
        folder_path = path.join(
            working_folder, 'event_%d' % folder_id, 
            'OSCAR_events', events_list[iev].split('/')[-1])
        bashCommand = "ln -s %s %s" % (
            path.abspath(events_list[iev]), folder_path)
        subprocess.Popen(bashCommand, stdout = subprocess.PIPE, shell=True)
        mixed_id = random.randint(0, len(events_list)-1)
        filename_mixed = events_list[mixed_id].split('/')[-1].split('.')[0]    
        mixed_event_id = filename_mixed.split('_')[-1]                         
        while (mixed_event_id == iev):                                         
            mixed_id = random.randint(0, len(events_list)-1)                   
            filename_mixed = events_list[mixed_id].split('/')[-1].split('.')[0]
            mixed_event_id = filename_mixed.split('_')[-1]                     
        folder_path = path.join(                                               
                    working_folder, 'event_%d' % folder_id,                            
                    'OSCAR_events', 'mixed_event_%s.dat' % event_id)                   
        bashCommand = "ln -s %s %s" % (                                        
                    path.abspath(events_list[mixed_id]), folder_path)                  
        subprocess.Popen(bashCommand, stdout = subprocess.PIPE, shell=True)    


def generate_event_folder_iSS(cluster_name, working_folder, event_id):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    mkdir(path.join(event_folder, 'hydro_events'))
    generate_script_iSS(cluster_name, event_folder)
    shutil.copytree('codes/iSS', 
                    path.join(path.abspath(event_folder), 'iSS'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/iSS_code/iSS_tables'),
        path.join(path.abspath(event_folder), "iSS/iSS_tables")), shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/iSS_code/iSS.e'),
        path.join(path.abspath(event_folder), "iSS/iSS.e")), shell=True)
    shutil.copytree('codes/osc2u', 
                    path.join(path.abspath(event_folder), 'osc2u'))
    shutil.copytree('codes/urqmd', 
                    path.join(path.abspath(event_folder), 'urqmd'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/urqmd_code/urqmd/urqmd.e'),
        path.join(path.abspath(event_folder), "urqmd/urqmd.e")), shell=True)
    shutil.copytree('codes/hadronic_afterburner_toolkit', 
                    path.join(path.abspath(event_folder), 
                    'hadronic_afterburner_toolkit'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.abspath('codes/hadronic_afterburner_toolkit_code/EOS'),
        path.join(path.abspath(event_folder),
                  "hadronic_afterburner_toolkit/EOS")), shell=True)

def generate_event_folder_iS(cluster_name, working_folder, event_id):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    mkdir(path.join(event_folder, 'hydro_events'))
    generate_script_iS(cluster_name, event_folder)
    shutil.copytree('codes/iS',
                    path.join(path.abspath(event_folder), 'iS'))

def copy_hydro_events(number_of_cores, input_folder, working_folder):
    events_list = glob('%s/surface*.dat' % input_folder)
    for iev in range(len(events_list)):
        event_id = events_list[iev].split('/')[-1].split('_')[-1].split('.')[0]
        folder_id = iev % number_of_cores
        working_path = path.join(working_folder, 'event_%d' % folder_id,
                                 'hydro_events')
        folder_path = path.join(working_path, events_list[iev].split('/')[-1])
        bashCommand = "ln -s %s %s" % (
            path.abspath(events_list[iev]), folder_path)
        subprocess.Popen(bashCommand, stdout = subprocess.PIPE, shell=True)
        shutil.copy(path.join(input_folder, 
                              'music_input_event_%s' % event_id), 
                    working_path)

def copy_job_scripts(working_folder):
    shutil.copy("job_MPI_wrapper.py", working_folder)
    shutil.copy("submit_MPI_job_for_all.pbs", working_folder)
    shutil.copy("run_job.sh", working_folder)

def print_mode_cheat_sheet():
    print("Here is a cheat sheet for mode option:")
    print("mode -1: run iS + resonance decay")
    print("mode 0: run iSS + osc2u + UrQMD from hydro hypersurface")
    print("mode 1: run UrQMD with OSCAR events")
    print("mode 2: calculate HBT correlation with OSCAR events")
    print("mode 3: calculate HBT correlation with UrQMD events")
    print("mode 4: collect spectra and flow observables from UrQMD events")
    print("mode 5: run JAM with OSCAR events")
    print("mode 6: collect spectra and vn with JAM events")
    print("mode 7: calculate HBT correlation with JAM events")
    print("mode 8: collect particle yield distribution with UrQMD events")
    print("mode 9: collect particle yield distribution with OSCAR events")

if __name__ == "__main__":
    try:
        from_folder = str(sys.argv[1])
        folder_name = str(sys.argv[2])
        cluster_name = str(sys.argv[3])
        ncore = int(sys.argv[4])
        mode = int(sys.argv[5])
    except IndexError:
        print("Usage:")
        print("  %s input_folder working_folder cluster_name num_of_cores mode"
              % str(sys.argv[0]))
        print("")
        print_mode_cheat_sheet()
        exit(0)

    if mode == 0:   # run iSS + osc2u + UrQMD from hydro hypersurface
        for icore in range(ncore):
            generate_event_folder_iSS(cluster_name, folder_name, icore)
        copy_hydro_events(ncore, from_folder, folder_name)
        copy_job_scripts(folder_name)
    elif mode == -1:   # run iS + resonance decay
        for icore in range(ncore):
            generate_event_folder_iS(cluster_name, folder_name, icore)
        copy_hydro_events(ncore, from_folder, folder_name)
    elif mode == 1:   # run UrQMD with OSCAR events
        for icore in range(ncore):
            generate_event_folder(cluster_name, folder_name, icore)
        copy_OSCAR_events(ncore, from_folder, folder_name)
    elif mode == 2:   # calculate HBT correlation with OSCAR events
        for icore in range(ncore):
            generate_event_folder_UrQMD(cluster_name, folder_name, icore, mode)
        copy_OSCAR_events(ncore, from_folder, folder_name)
    elif mode == 3:   # calculate HBT correlation with UrQMD events
        for icore in range(ncore):
            generate_event_folder_UrQMD(cluster_name, folder_name, icore, mode)
        copy_UrQMD_events(ncore, from_folder, folder_name)
        copy_job_scripts(folder_name)
    elif mode == 4:   # collect spectra and flow observables from UrQMD events
        for icore in range(ncore):
            generate_event_folder_UrQMD(cluster_name, folder_name, icore, mode)
        copy_UrQMD_events(ncore, from_folder, folder_name)
        copy_job_scripts(folder_name)
    elif mode == 5:   # run JAM with OSCAR events
        for icore in range(ncore):
            generate_event_folder_JAM(cluster_name, folder_name, icore, mode)
        copy_OSCAR_events(ncore, from_folder, folder_name)
    elif mode == 6:   # collect spectra and vn with JAM events
        for icore in range(ncore):
            generate_event_folder_JAM(cluster_name, folder_name, icore, mode)
        copy_JAM_events(ncore, from_folder, folder_name)
    elif mode == 7:   # calculate HBT correlation with JAM events
        for icore in range(ncore):
            generate_event_folder_JAM(cluster_name, folder_name, icore, mode)
        copy_JAM_events(ncore, from_folder, folder_name)
    elif mode == 8:  # collect particle yield distribution with UrQMD events
        for icore in range(ncore):
            generate_event_folder_UrQMD(cluster_name, folder_name, icore, mode)
        copy_UrQMD_events(ncore, from_folder, folder_name)
    elif mode == 9:  # collect particle yield distribution with OSCAR events
        for icore in range(ncore):
            generate_event_folder_UrQMD(cluster_name, folder_name, icore, mode)
        copy_OSCAR_events(ncore, from_folder, folder_name)
    elif mode == 10:   # calculate balance function correlation with UrQMD events
        for icore in range(ncore):
            generate_event_folder_UrQMD(cluster_name, folder_name, icore, mode)
        copy_UrQMD_events(ncore, from_folder, folder_name)
        copy_job_scripts(folder_name)


