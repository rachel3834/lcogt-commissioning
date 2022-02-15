import argparse
import datetime
import logging
import tempfile
import threading
from ctypes import *
import numpy as np
import time
import astropy.io.fits as fits
import vxi11
import requests
import time
from io import BytesIO
_logger = logging.getLogger(__name__)

requests_logger = logging.getLogger('connectionpool')
requests_logger.setLevel(logging.ERROR)

class LCOLab:

    def __init__(self):
        self.ins = vxi11.Instrument('172.16.4.6')

        _logger.debug ("VXI11 interface: %s" % (self.ins.ask("*IDN?")))

    def expose(self, exptime, voltage = None, block=True, overhead=0):
        ''' Expose via constant illumination.
        Sets the function generator into pulse mode'''
        _logger.debug ("Lab exposing for % 5.2f s" % (exptime))
        self.ins.write ("puls:mode TRIG")
        self.ins.write (f"burst:ncycles 1")
        self.ins.write ("puls:state ON")
        self.ins.write ("puls:per %fs" % (exptime+overhead+1))
        self.ins.write ("puls:widt %fs" % (exptime+overhead))
        self.ins.write ("PULSe:DELay 0s")
        self.ins.write ("burst:DELay 0s")

        # self.ins.write (f"PULSe:DCYCs 100" )
        if (voltage is not None) and (voltage >=0.) and (voltage <= 5.):
            _logger.debug (f"Setting LED voltage to {voltage}")
            self.ins.write(f"voltage:level:imm:high {voltage} V")

        self.ins.trigger()
        if block:
            _logger.info("Blocking during exposure time")
            time.sleep (exptime)
        _logger.debug ("Done exposing")


    def expose_burst (self, exptime,  frequency=100, ncycles = 10, voltage=None, block=True, overhead=5):
        _logger.info (f"Lab burst exposing for {exptime}, led {voltage}, overhead {overhead}")
        time.sleep(overhead)
        _logger.info ("Done sleeping, firing up the LED")
        if ncycles == 0:
            return
        self.ins.write ("burst:state ON")
        self.ins.write ("burst:mode TRIG")
        self.ins.write (f"burst:ncycles {ncycles}")
        self.ins.write (f"freq:fixed {frequency}Hz")
        self.ins.write (f"PULSe:DCYC 50.0" )
        self.ins.write (f"burst:DELay {overhead}s")
        self.ins.write (f"pulse:DELay {overhead}s")

        if (voltage is not None) and (voltage >=0.) and (voltage <= 5.):
            _logger.debug (f"Setting LED voltage to {voltage}")
            self.ins.write(f"voltage:level:imm:high {voltage} V")

        self.ins.trigger()
        if block:
            _logger.info("Blocking during exposure time")
            time.sleep (exptime)
            _logger.debug ("Done exposing")



    def close(self):
        self.ins.close()


class restcam:


    def __init__(self, ipaddr):
        self.ipaddr = ipaddr
        self.bin_x = 1
        self.bin_y = 1

    """ Relase camera and close sdk """
    def close(self):
        pass


    def __str__(self):
        w,h = self.getPixelSize()
        return f" {str(self.id.value.decode('UTF-8'))} Pixels: {w} x {h} pixels"

    def setReadMode (self, readmodenumber):
        pass


    def sendSettings (self, payload):
        r = requests.get(f"http://{self.ipaddr}/api/ImagerSetSettings.cgi", data = payload)
        return r

    def setBinning(self, wbin, hbin):
        self.bin_x = wbin
        self.bin_y = hbin
        payload = {'BinX' : self.bin_x,
                   'BinY' : self.bin_y}
        r = self.sendSettings(payload)
        _logger.info(f"Binning returned {r}")


    def setGain(self, gain):
        """ Set camera gain """
        self.gain = gain

    def getGain (self):
        return 1.0


    def setTemp (self, temperature):
        self.setpoint = (temperature)




    def getTemp (self):

        return 0.0



    def setROI(self,x,y):

        pass



    def getChipDimensions(self):
        return self.chipw.value, self.chiph.value

    def getPixelSize(self):
        return self.w.value, self.h.value

    def getStatus (self):
        try:
            status = requests.get (f"http://{self.ipaddr}/api/ImagerState.cgi").json()
        except:
            logging.error(f"Error while asking for SBIG status : {status}")
        return status

    def getframe(self, exptime, filename, args=None):

        _logger.info (f"Starting exposure {exptime} seconds")

        frametype = 1
        if ('b00' in filename) or ('d00') in filename:
            frametype=2

        params = {
            'Duration': exptime,
            'FrameType' : frametype,
        }
        r = requests.get (f"http://{self.ipaddr}/api/ImagerStartExposure.cgi", params = params  )
        _logger.debug(f"Started exposure {r}")
        time.sleep (0.5)

        while ( self.getStatus() > 0):
            time.sleep (0.5)
        start_readout = datetime.datetime.utcnow()
        _logger.info ("Downloading image from camera")
        fitsdata = requests.get (f"http://{self.ipaddr}/api/Imager.FIT", params ={}  )

        _logger.info (f"writing to fits file {filename}")
        with open(filename, "wb") as myfile:
            myfile.write (fitsdata.content)
            myfile.close()

        hdul = fits.open (filename)

        object = None
        if args is not None:
            object = f"led {args.ledvoltage} nburst {args.nburstcycles}"

        prihdr = hdul[0].header
        prihdr['OBJECT'] = object
        prihdr['EXPTIME'] = exptime
        prihdr['FILTER'] = 'None'
        prihdr['AIRMASS'] = 1.0
        prihdr['DATE-OBS'] = start_readout.strftime('%Y-%m-%dT%H:%M:%S')
        #prihdr['GAIN'] = self.gain.value
        prihdr['CCDSUM'] = f"[{self.bin_x} {self.bin_y}]"
        prihdr['READMODE'] = 'default'
        prihdr['CCDTEMP'] = prihdr['CCD-TEMP']

        hdul.writeto (filename,  overwrite=True )
        end_fitswrite = datetime.datetime.utcnow()
        _logger.info("done taking image")




def main():
    args = parseCommandLine()


    if args.testled is not None:
        lab = LCOLab()
        print ("LED ON")
        lab.expose(exptime = args.testled, overhead = 1, block=True, voltage=args.ledvoltage)
        print ("LED OFF")
        exit (0)

    if args.testledburst:
        lab = LCOLab()
        print ("LED ON")
        lab.expose_burst(exptime = args.testledburst, overhead = 1, block=True, voltage=args.ledvoltage)
        print ("LED OFF")
        exit (0)

    qhyccd = restcam("10.6.11.64:8080")


    if args.settemp is not None:
        qhyccd.setTemp(args.settemp)
        qhyccd.close()
        exit(0)


    if args.gettemp:
        qhyccd.getTemp()
        qhyccd.close()
        exit(0)

    qhyccd.setBinning (args.binning, args.binning)

    suffix = None
    if args.bias:
        suffix = 'b00'
    if args.dark:
        suffix = 'd00'
    if args.flat:
        suffix = 'f00'
        lab = LCOLab()


    for exptime in args.exptime:
        _logger.info (f"taking exposures for exptime {exptime}")
        for ii in range (args.expcnt):
            imagename=f"{args.outputpath}/restcam-{datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S')}.{suffix}.fits"
            if args.flat and lab is not None:
                if args.nburstcycles is None:
                    # This is a conventional exposure where we ensure the LED is on befor we open the shutter and stays on until shutter closes.
                    _logger.info ("Starting conventional shutter-defined exposure")
                    lab.expose(exptime = exptime, overhead = 2, block=False, voltage=args.ledvoltage)
                else:
                    # Here we open the shutter, and then turn the LED on for a determined amount of time. it takes a few seconds from requesting an exposure
                    # until the shutter actually opens. Hence we are putting the LED con command into a background thread that starts its working day with sleeping.

                    _logger.info (f"Starting frequencey generator defined exposure for {args.nburstcycles} cycles.")
                    th =threading.Thread ( target=lab.expose_burst, kwargs={'exptime':exptime, 'ncycles':args.nburstcycles, 'overhead':7, 'voltage':args.ledvoltage, 'block':False})
                    th.start()

            qhyccd.getframe(exptime, imagename, args=args)

            if args.flat:
                time.sleep (1)


    qhyccd.close()
    if args.flat:
        lab.close()
    exit(0)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='QHYCCD command line control')

    actions = parser.add_mutually_exclusive_group()
    actions.add_argument ("--bias", action = "store_true", help="Take bias exposures. Set number via --expcnt" )
    actions.add_argument ("--dark", action = "store_true", help="Take dark exposures. Set number via --expcnt and darktime via --exptime" )
    actions.add_argument ("--flat", action = "store_true", help="Take flat exposures. Set number via --expcnt and darktime via --exptime" )
    actions.add_argument ("--settemp", type=float, help="Set CCD target temperature" )
    actions.add_argument ("--gettemp", action = "store_true",  help="get CCD target temperature" )
    actions.add_argument ("--testled", type=float,  help="testled" )
    actions.add_argument ("--testledburst", type=float,  help="testled in burst mode" )

    actions.add_argument ("--chamberpump", type = bool,  help="cycle detector chaber decissitant" )


    parser.add_argument ("--ledvoltage", type=float, default = None)
    parser.add_argument('--expcnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, nargs='*', default=[0,])
    parser.add_argument('--gain', type=int, default=5)
    parser.add_argument('--binning', type=int, default=1)
    parser.add_argument('--readmode', type=int, default=0)
    parser.add_argument('--nburstcycles', type=int, default = None)



    parser.add_argument('--outputpath', type=str, default="data", help="outputpath")
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    return args





if __name__ == '__main__':
    main()
