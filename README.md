# lcogt-commissioning
Software for use in astronomical commissioning of instruments on the LCOGT network

Updated by Daniel Harbeck:
==
The original scripts provided by Rachel Street were tuned for the old 3D cube fits format
 of the LCO Sinistro cameras. Since the original commissioning of the Sinsitro cameras, the storage
 format was revised and now uses a MEF format.For continued commissioning and reverification of the
 LCO cameras, these tools are being rewritten by Daniel Harbeck to be compatible with the new format. 

 Note: The scripts explained below are updated; scripts not mentioned here have not been vetted yet and 
 is vestigial code 

noisegainmef.py
===

Measure read noise [e-] and gain [e-/ADU]  based on a pair of raw biases images and equally exposed
flat field images. 

Syntax is:
<pre>
usage: noisegainrawmef.py [-h] [--imagepath OPT_IMAGEPATH]
                          [--log_level {DEBUG,INFO,WARN}]
                          fitsfile fitsfile fitsfile fitsfile

General purpose noise and gain measurement from a set of two flat fields and
two biases.

positional arguments:
  fitsfile              four fits files: bias_1 bias_2 flat_1 flat_2

optional arguments:
  -h, --help            show this help message and exit
  --imagepath OPT_IMAGEPATH
                        pathname to prepend to fits file names.
  --log_level {DEBUG,INFO,WARN}
                        Set the debug level

</pre>

The fitsfiles are to be in the order: bias1, bias2, flat1, flat2

Example:
<pre>
python noisegainrawmef.py --imagepath /nfs/archive/engineering/cpt/fl14/20170912/raw/ cpt1m013-fl14-20170912-0002-b00.fits.fz cpt1m013-fl14-20170912-0003-b00.fits.fz   cpt1m013-fl14-20170912-0030-f00.fits.fz cpt1m013-fl14-20170912-0031-f00.fits.fz
</pre>


crosstalk.py
====
Measure the crosstalk between amplifiers in a MEF file. 

Updated command line parsing to use argparse. Updated fits reading to properly handle .fz
 compressed fits files. 

Usage example:
<pre>

python crosstalk.py --quadrant 2 --linear --imagepath /nfs/archive/engineering/cpt/fl14/20170913/raw    cpt1m013-fl14-20170913-0106-x00.fits.fz cpt1m013-fl14-20170913-0107-x00.fits.fz cpt1m013-fl14-20170913-0108-x00.fits.fz
</pre>

Paramters:

--quadrant [N] parameter defines which extension number (start counting at 0 ) is considered the contaminating quadrant.
  When loading a raw image in ds9 with the -mosaicimage iraf option, this is where the quadrants are:
<pre>  
       4|3
       -+- 
       1|2
</pre>
  the input parameter needs those quadrant ids -1 (should change code one day though!)

--linear will determine linear crosstalk relation

--imagepath [path]  points to the directory where images located and will be prepended to the
 following individual image names. Save a bit of repetitive typing of the same path.

The following fits images wil be analysed for cross talk. While the software was written with the 
LCO Sinistro cameras in mind, the current version o this script should work with any camera as 
long as the output of amplifiers is stored in equally dimensioned fits extensions. 