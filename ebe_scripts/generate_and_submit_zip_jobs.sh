#!/usr/bin/env bash

working_folder=$1

(
    ./zip_urqmd_events.py $working_folder
    cp ./codes/hadronic_afterburner_toolkit_code/convert_to_binary.e $working_folder
    cd $working_folder
    qsub *.pbs
)
