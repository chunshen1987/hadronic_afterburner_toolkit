#!/usr/bin/env bash

fromFolder=$1
toFolder=$2

mkdir $toFolder
cp $fromFolder/*.f13 $toFolder
(
cd $toFolder
for ii in `ls`
do
   mv $ii particle_list_`echo $ii | cut -f 2 -d _ | cut -f 1 -d .`.dat
done
)
