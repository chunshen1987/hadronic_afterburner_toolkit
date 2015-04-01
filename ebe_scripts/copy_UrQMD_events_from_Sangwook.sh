#!/usr/bin/env bash

fromFolder=$1
toFolder=$2

mkdir $toFolder
(
cd $toFolder
for ii in `ls $fromFolder | grep "f13"`
do
   ln -s $fromFolder/$ii particle_list_`echo ${ii##*/} | cut -f 2 -d _ | cut -f 1 -d .`.dat
done
)
