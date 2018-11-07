#!/usr/bin/env bash

usage="./collect_events.sh fromFolder toFolder"

fromFolder=$1
toFolder=$2

if [ -z "$fromFolder" ]
then
    echo $usage
    exit 1
fi
if [ -z "$toFolder" ]
then
    echo $usage
    exit 1
fi

echo "collecting events from " $fromFolder " to " $toFolder

folderName=`echo $fromFolder | sed 's/arena_//'`
mkdir $toFolder/$folderName

total_eventNum=0
collected_eventNum=0
for ijob in `ls --color=none $fromFolder`;
do 
    eventsPath=$fromFolder/$ijob/spvn_results
    for iev in `ls --color=none $eventsPath`
    do 
        if [ -a $eventsPath/$iev/particle_9999_vndata_eta_-0.5_0.5.dat ]; then
            mv $eventsPath/$iev $toFolder/$folderName
            ((collected_eventNum++))
        fi
        ((total_eventNum++))
    done
done

echo "Collected events number: " $collected_eventNum " out of " $total_eventNum

./combine_results_into_hdf5.py $toFolder/$folderName
mv $folderName.h5 $toFolder/
rm -fr $toFolder/$folderName
