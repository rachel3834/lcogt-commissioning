#!/bin/bash

# Crawls thourgh the LCO archive and attempts to fit gain / readnosie from automatically identified pairs of sky flats
# and daytime biases.
# Edit according to your need below.
# Reqruiements: gnu parallel needs to be installed.



noisegaindatabase="noisegain.sqlite"   # Output database file to store results.
webpageoutputdir="gainhistory"         # existing local subdiretory into which an overview html page will be rendered.


base=/archive/engineering              # Location of archive mount
cameras="fa?? fs?? kb??"               # Camera types to crawl.
sites="elp cpt lsc ogg coj tfn"        # Which sites to crawl
dates="201905??"                       # date range selector.

inputselection="*-[bf]00.fits.fz"      # type of images to serach for. need flats and biases here.

NCPU=30                                # How many CPUs to employ

for site in $sites; do
 for camera in $cameras; do

  sitecameras=`find ${base}/${site}  -maxdepth 1 -type d -wholename "*/$camera"`

  for sitecamera in $sitecameras; do

   directories=`find "${sitecamera}" -maxdepth 1 -type d  -wholename "*/${dates}" `

   for day in $directories; do

     searchpath=${day}/raw/${inputselection}
     echo "Searchpath is $searchpath"
     sem  -j $NCPU python lcocommissioning/noisegainrawmef.py --noreprocessing --loglevel INFO --sortby filterlevel --database ${noisegaindatabase} $searchpath

   done

  done

 done
done

sem --wait

analysgainhistory --outputdir ${webpageoutputdir} --database ${noisegaindatabase}
