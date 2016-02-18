#!/usr/bin/env bash

fromFolder=$1
toFolder=$2

target_surface_filename="surface_eps_0.1828.dat"

mkdir $toFolder
(
  cd $toFolder
  for iev in `ls --color=none $fromFolder | grep "event-"`
  do
     event_id=`echo $iev | sed 's/event-//'`
     ln -s $fromFolder/$iev/$target_surface_filename surface_event_$event_id.dat
     cp $fromFolder/$iev/music_input music_input_event_$event_id
  done
)
