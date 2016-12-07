#!/usr/bin/env python

import sys
from os import path

def generate_job_script(folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    walltime = '10:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    script.write(
"""#!/usr/bin/env bash
#PBS -N zip_urqmd_events
#PBS -l nodes=1:ppn=1
#PBS -l walltime=%s
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -A cqn-654-ad
#PBS -q sw
#PBS -d %s

for ii in `ls *.dat`
do
    ./convert_to_binary.e $ii
    rm -fr $ii
done

""" % (walltime, working_folder))
    script.close()

if __name__ == "__main__":
    try:
        folder_name = str(sys.argv[1])
    except IndexError:
        print "%s folder_name" % str(sys.argv[0])
        exit(0)
    generate_job_script(folder_name)

