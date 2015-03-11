#!/usr/bin/env python

import sys
from os import path, mkdir
import shutil
from glob import glob

def generate_script(folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    walltime = '1:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
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
#PBS -m bea
#PBS -M chunshen1987@gmail.com
#PBS -d %s

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

""" % (working_folder.split('/')[-1], walltime, working_folder))
    script.close()

def generate_script_HBT(folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    walltime = '2:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    script.write(
"""#!/usr/bin/env bash
#PBS -N HBT_%s
#PBS -l nodes=1:ppn=1
#PBS -l walltime=%s
#PBS -S /bin/bash
#PBS -A cqn-654-ad
#PBS -e test.err
#PBS -o test.log
#PBS -q sw
#PBS -m bea
#PBS -M chunshen1987@gmail.com
#PBS -d %s

mkdir HBT_results
for iev in `ls UrQMD_events`
do
    cd HBTcorrelation_MCafterburner
    rm -fr results
    mkdir results
    mv ../UrQMD_events/$iev results/particle_list.dat
    ./HBT_afterburner.e > output.log
    mv results/particle_list.dat ../UrQMD_events/$iev
    mv results ../HBT_results/event_`echo $iev | cut -f 3 -d _ | cut -f 1 -d .`
    cd ..
done

""" % (working_folder.split('/')[-1], walltime, working_folder))
    script.close()

def generate_script_HBT_with_OSCAR(folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    script.write(
"""#!/usr/bin/env bash
#PBS -N HBT_%s
#PBS -l nodes=1:ppn=1:shai-hulud
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -m bea
#PBS -M chunshen1987@gmail.com
#PBS -d %s

mkdir HBT_results
for iev in `ls OSCAR_events`
do
    cd HBTcorrelation_MCafterburner
    rm -fr results
    mkdir results
    mv ../OSCAR_events/$iev results/OSCAR.DAT
    ./HBT_afterburner.e mode=0 > output.log
    mv results/OSCAR.DAT ../OSCAR_events/$iev
    mv results ../HBT_results/event_`echo $iev | cut -f 2 -d _ | cut -f 1 -d .`
    cd ..
done

""" % (working_folder.split('/')[-1], working_folder))
    script.close()

def copy_UrQMD_events(number_of_cores, working_folder):
    events_list = glob('UrQMD_results/particle_list_*.dat')
    for iev in range(len(events_list)):
        folder_id = iev % number_of_cores
        folder_path = path.join(working_folder, 'event_%d' % folder_id)
        shutil.move(events_list[iev], 
                    path.join(folder_path, 'UrQMD_events', 
                              events_list[iev].split('/')[-1]))


def generate_event_folder_UrQMD(working_folder, event_id):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)

    # calucalte HBT correlation with OSCAR outputs
    #mkdir(path.join(event_folder, 'OSCAR_events'))
    #generate_script_HBT_with_OSCAR(event_folder)

    # calculate HBT correlation with UrQMD outputs
    mkdir(path.join(event_folder, 'UrQMD_events'))
    generate_script_HBT(event_folder)

    shutil.copytree('codes/HBTcorrelation_MCafterburner', 
                    path.join(path.abspath(event_folder), 
                    'HBTcorrelation_MCafterburner'))


def generate_event_folder(working_folder, event_id):
    event_folder = path.join(working_folder, 'event_%d' % event_id)
    mkdir(event_folder)
    mkdir(path.join(event_folder, 'OSCAR_events'))
    generate_script(event_folder)
    shutil.copytree('codes/osc2u', 
                    path.join(path.abspath(event_folder), 'osc2u'))
    shutil.copytree('codes/urqmd', 
                    path.join(path.abspath(event_folder), 'urqmd'))


def copy_OSCAR_events(number_of_cores, working_folder):
    events_list = glob('OSCAR_events/*.dat')
    for iev in range(len(events_list)):
        folder_id = iev % number_of_cores
        folder_path = path.join(working_folder, 'event_%d' % folder_id)
        shutil.move(events_list[iev], 
            path.join(folder_path, 'OSCAR_events', 
                      events_list[iev].split('/')[-1]))

if __name__ == "__main__":
    folder_name = str(sys.argv[1])
    ncore = int(sys.argv[2])
    mode = int(sys.argv[3])
    if mode == 1:
        for icore in range(ncore):
            generate_event_folder(folder_name, icore)
        copy_OSCAR_events(ncore, folder_name)
    elif mode == 2:
        for icore in range(ncore):
            generate_event_folder_UrQMD(folder_name, icore)
        # calculate HBT correlation from OSCAR files
        #copy_OSCAR_events(ncore, folder_name)
        # calculate HBT correlation from UrQMD files
        copy_UrQMD_events(ncore, folder_name)
