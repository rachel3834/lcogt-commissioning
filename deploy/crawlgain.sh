#!/bin/bash

# Crawls thourgh the LCO archive and attempts to fit gain / readnosie from automatically identified pairs of sky flats
# and daytime biases.
# Edit according to your need below.
# Reqruiements: gnu parallel needs to be installed.



noisegaindatabase="sqlite:///noisegain.sqlite"   # Output database file to store results.
webpageoutputdir="/home/dharbeck/public_html/gainhistory"         # existing local subdiretory into which an overview html page will be rendered.
reprocessing="--noreprocessing"

base=/archive/engineering              # Location of archive mount

NDAYS=3

time crawlnoisegain --ndays $NDAYS $reprocessing  --cameratype "fa??" --readmode "full_frame" --loglevel INFO --database ${noisegaindatabase}
time crawlnoisegain --ndays $NDAYS $reprocessing  --cameratype "fa??" --readmode "central_2k_2x2" --loglevel INFO --database ${noisegaindatabase}
time crawlnoisegain --ndays $NDAYS $reprocessing --cameratype "fs??" --readmode "default" --loglevel INFO --database ${noisegaindatabase}
time crawlnoisegain --ndays $NDAYS $reprocessing --cameratype "kb??" --readmode "default" --loglevel INFO --database ${noisegaindatabase}

time analysegainhistory --outputdir ${webpageoutputdir} --database ${noisegaindatabase}
