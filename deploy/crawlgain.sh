#!/bin/bash

noisegaindatabase="noisegain.sqlite"
webpageoutputdir="gainhistory"

base=/archive/engineering
cameras="fa?? fl?? fs?? kb??"
sites="elp cpt lsc ogg coj tfn"
dates="2018????"

inputselection="*-[bf]00.fits.fz"

NCPU=2

for site in $sites; do
 for camera in $cameras; do

  sitecameras=`find ${base}/${site}  -maxdepth 1 -type d -wholename "*/$camera"`

  for sitecamera in $sitecameras; do

   directories=`find "${sitecamera}" -maxdepth 1 -type d  -wholename "*/${dates}" `

   for day in $directories; do

     searchpath=${day}/raw/${inputselection}
     echo "Searchpath is $searchpath"
     sem  -j $NCPU python noisegainrawmef.py --noreprocessing --loglevel INFO --sortby filterlevel --database ${noisegaindatabase} $searchpath

   done

  done

 done
done

sem --wait

python analysegainhistory.py --outputdir ${webpageoutputdir} --database ${noisegaindatabase}