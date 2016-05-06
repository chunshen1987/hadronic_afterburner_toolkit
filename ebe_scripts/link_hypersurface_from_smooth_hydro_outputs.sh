#!/usr/bin/env bash

fromFolder=$1
toFolder=$2

target_surface_filename="surface_eps_0.1.dat"

mkdir $toFolder
(
  cd $toFolder
  for iev in {0..29}
  do
     ln -s $fromFolder/$target_surface_filename surface_event_$iev.dat
     cp $fromFolder/music_input music_input_event_$iev
  done
)
