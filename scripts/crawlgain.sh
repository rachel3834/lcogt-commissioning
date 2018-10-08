#!/bin/bash

base=/archive/engineering/
camera=fl16
site=cpt

inputselection="*-[bf]00.fits.fz"


directories=`find ${base}/${site}/${camera} -type d -wholename "*/2018*/raw" `

for day in $directories
do
   echo $day
   searchpath=$day/$inputselection
   python noisegainrawmef.py --log_level DEBUG --sortby filterlevel $searchpath
done