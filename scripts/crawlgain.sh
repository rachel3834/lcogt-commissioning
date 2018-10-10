#!/bin/bash

base=/archive/engineering/
camera=fl04
site=lsc

inputselection="*-[bf]00.fits.fz"


directories=`find ${base}/${site}/${camera} -type d -wholename "*/201[67]*/raw" `

for day in $directories
do
   echo $day
   searchpath=$day/$inputselection
   python noisegainrawmef.py --log_level INFO --sortby filterlevel $searchpath
done