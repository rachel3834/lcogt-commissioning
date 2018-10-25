#!/bin/bash

base=/archive/engineering
cameras="fa?? fl?? fs?? kb??"
sites="elp cpt lsc ogg coj tfn"
dates="201809??"

inputselection="*-[bf]00.fits.fz"

NCPU=4

for site in $sites; do
 for camera in $cameras; do

  sitecameras=`find ${base}/${site}  -maxdepth 1 -type d -wholename "*/$camera"`

  for sitecamera in $sitecameras; do

   directories=`find "${sitecamera}" -maxdepth 1 -type d  -wholename "*/${dates}" `

   for day in $directories; do

     searchpath=${day}/raw/${inputselection}
     echo "Searchpath is $searchpath"
     sem  -j 2 python noisegainrawmef.py --log_level INFO --sortby filterlevel $searchpath

   done

  done

 done
done

sem --wait