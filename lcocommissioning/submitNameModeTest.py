import argparse
import json
import logging
import ephem
import math
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import datetime as dt
import  lcocommissioning.common.common as common



_logger = logging.getLogger(__name__)

goodXTalkTargets = ['auto', '91 Aqr', 'HD30562', '15 Sex', '30Psc', '51Hya']





def createRequestsForStar(context):
    exposuretime = context.exptime   # in minutes
    overheadperexposure = 7
    telescopeslew = 60
    start = context.start
    nexposure = int(context.expcnt)

    # create one block per quadrant
    end = start + dt.timedelta(seconds = nexposure * (exposuretime+overheadperexposure)) \
          + dt.timedelta(seconds=telescopeslew)
    start = str(start).replace(' ','T')
    end = str(end).replace(' ','T')

    offsetPointing = SkyCoord(context.radec.ra + Angle(context.offsetRA, unit=u.arcsec),
                             context.radec.dec + Angle(context.offsetDec, unit=u.arcsec))

    pointing =  {
        "type": "SP",
        "name": "Named Mode Test {}".format (context.name),
        "coord_type": "RD",
        "coord_sys": "ICRS",
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "pro_mot_ra": "0",
        "pro_mot_dec": "0",
        "parallax": "0.0000000",
        "ra" :  "%10f" % offsetPointing.ra.degree,
        "dec" : "%7f" % offsetPointing.dec.degree,
    }

    print("Block for {} from {} to {}".format ( pointing['name'], str(start), str(end)))

    block = {"molecules":[],
             'start': start,
             'end': end,
             'site': context.site,
             'observatory': context.dome,
             'telescope': context.telescope,
             'instrument_class': '1m0-SciCam-Sinistro'.upper(),
             'priority': 32,
             "is_too": False,
             "max_airmass": "3.0000000",
             "min_lunar_dist": "30.0000000",
             }

    for readoutmode in ("lco2_500kHz_binned_window",): # "","lco2_500kHz_DRH",):

        molecule = {
            "filter": context.filter,
            "exposure_time": exposuretime,
            "readout_mode": readoutmode,
            "pointing": pointing,
            "defocus": "0.0000000",
            "ag_mode": "OPT",
            "tracking_num": "",
            "request_num": "",
            "tag_id": "TST",
            "user_id": "daniel_harbeck",
            "prop_id": "LCOEngineering",
            "group": "test_group",
            "exposure_count": nexposure,
            "bin_x": 1,
            "bin_y": 1,
            "inst_name": context.instrument,
            "priority": 31,
            "type": "EXPOSE",
            "ag_filter": "",
            "ag_exp_time": "10.0000000"

        }



        block['molecules'].append (molecule)

        _logger.debug (json.dumps(block,indent=4))
        common.send_to_lake(block, context.opt_confirmed)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Submiot a named mode observation request to POND')

    parser.add_argument('--name', default='auto', type=str,
                        help='Name of target.. Will be resolved via simbad. If resolve failes, program will exit.')
    parser.add_argument('--defocus', type=float, default=0.0, help="Amount to defocus star.")
    parser.add_argument('--site',  choices=common.lco_1meter_sites, required=True,
                        help="To which site to submit")
    parser.add_argument('--dome',  choices=['doma', 'domb', 'domc'], required=True, help="To which enclosure to submit")
    parser.add_argument('--telescope', default='1m0a',)
    parser.add_argument('--filter', default='rp')
    parser.add_argument('--exp-cnt',  type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime',  type=float, default=20)
    parser.add_argument('--instrument', default=common.lco_sinistro1m_cameras,
                        help="To which instrument to submit")
    parser.add_argument('--start', default=None,
                        help="When to start observation. If not given, defaults to \"NOW\"")
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")

    # Per default, do not be on chip gap!
    parser.add_argument('--offsetRA', default=60, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=60, help="Extra pointing offset to apply Dec")

    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, block will be submitted.')

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    if args.start is None:
        args.start = dt.datetime.utcnow()
    else:
        try:
            args.start = dt.datetime.strptime(args.start, "%Y%m%d %H:%M")
        except ValueError:
            _logger.error("Invalidt start time argument: ", args.start)
            exit(1)

    if ('auto' in args.name):
        args.name = common.get_auto_target(goodXTalkTargets, args.site, args.start)
        if args.name is None:
            exit (1)

    try:
        _logger.debug("Resolving target name")
        args.radec = SkyCoord.from_name(args.name)
    except:
        print("Resolving target name failed, giving up")
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.name, args.radec.ra, args.radec.dec))
    return args


def main():
    args = parseCommandLine()
    createRequestsForStar(args)

if __name__ == '__main__':
   main()
