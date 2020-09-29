import argparse
import copy
import json
import logging
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import datetime as dt
import lcocommissioning.common.common as common

_logger = logging.getLogger(__name__)

goodFloydsFluxStandards = ['auto', 'HZ 43', 'GD 71', 'BD+284211', 'HZ 44', 'L745-46A', 'Feige 110', 'EGGR274']


def createRequestsForStar(context):
    exposuretime = context.exptime
    overheadperexposure = 60  # readout plus fixed overhead
    telescopeslew = 120 + 90 + 60  # telscope slew, acquire exposure + process, config change time
    start = context.start
    nexposure = int(context.expcnt)

    # create one block per quadrant
    end = start + \
          dt.timedelta(seconds=telescopeslew) + \
          dt.timedelta(seconds=nexposure * (exposuretime + overheadperexposure)) + \
          dt.timedelta(seconds=2 * (80 + overheadperexposure) + (80 + overheadperexposure))  # 2x flat, 1x arc

    start = str(start).replace(' ', 'T')
    end = str(end).replace(' ', 'T')

    offsetPointing = SkyCoord(context.radec.ra + Angle(context.offsetRA, unit=u.arcsec),
                              context.radec.dec + Angle(context.offsetDec, unit=u.arcsec))

    data = {
        'name': "Floyds test {}".format(context.targetname),
        'proposal': 'LCOEngineering',
        'site': context.site,
        'enclosure': 'clma',
        'telescope': '2m0a',
        'start': start,
        'end': end,
        'request': {
            'acceptability_threshold': 100,
            'configurations': [
                {'type': 'SPECTRUM',
                 'instrument_type': '2M0-FLOYDS-SCICAM',
                 'target': {
                     'type': 'ICRS',
                     'name': context.targetname,
                     'ra': offsetPointing.ra.degree,
                     'dec': offsetPointing.dec.degree
                 },
                 'acquisition_config': {
                     'mode': 'WCS'
                 },
                 'guiding_config': {
                     'mode': 'ON',
                     'optional': False
                 },
                 'constraints': {},
                 'instrument_configs': [{
                     'exposure_time': context.exptime,
                     'exposure_count': int(context.expcnt),
                     'mode': 'default',
                     'rotator_mode': 'VFLOAT',
                     'optical_elements': {
                         'slit': 'slit_1.2as'
                     }
                 }]
                 }
            ]
        }
    }


    # TODO: move to a better location in code!
    # agname = 'kb42' if 'ogg' in context.site else 'kb38'


    _logger.info(json.dumps(data, indent=4))
    common.submit_observation(data, context.opt_confirmed)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Submit an engineering Floyds observation.')

    parser.add_argument('--targetname', default='auto', type=str,
                        help='Name of star for Floyds test observation. if none is given, or auto, Software will try to'
                             ' find a cataloged flux stadnard star. Name must resolve via Simbad; if resolve failes,'
                             ' program will exit.')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--site', choices=common.lco_2meter_sites, required=True, help="To which site to submit")

    parser.add_argument('--exp-cnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, default=150)
    parser.add_argument('--slit', type=str, default="slit_1.2as", choices=['slit_1.2as', 'slit_2.0as', 'slit_6.0as'])

    parser.add_argument('--start', default=None,
                        help="When to start Floyds observation. If not given, defaults to \"NOW\"")
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")

    # Per default, do not be on chip gap!
    parser.add_argument('--offsetRA', default=0, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=0, help="Extra pointing offset to apply Dec")

    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, block will be submitted.')

    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    args.instrument = 'floyds01' if "ogg" in args.site else 'floyds02'
    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    if args.start is None:
        args.start = dt.datetime.utcnow()
    else:
        try:
            args.start = dt.datetime.strptime(args.start, "%Y%m%d %H:%M")
        except ValueError:
            _logger.error("Invalid start time argument: ", args.start)
            exit(1)

    if ('auto' in args.targetname):
        args.targetname = common.get_auto_target(goodFloydsFluxStandards, args.site, args.start)
        if args.targetname is None:
            exit(1)

    try:
        _logger.debug("Resolving target name")
        args.radec = SkyCoord.from_name(args.targetname)
    except:
        print("Resolving target name failed, giving up")
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.targetname, args.radec.ra, args.radec.dec))
    return args


def main():
    args = parseCommandLine()
    createRequestsForStar(args)
    print("Need to update this program for direct submission. Quitting.")
    exit(1)


if __name__ == '__main__':
    main()
