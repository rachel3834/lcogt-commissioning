import argparse
import copy
import json
import logging
import ephem
import math
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import datetime as dt
import requests

LAKE_URL = 'http://lake.lco.gtn'

_logger = logging.getLogger(__name__)

quadrantOffsets = {0: [-450, 450],
                   1: [450, 450],
                   2: [450, -450],
                   3: [-450, -450]}

_site_lonlat = {}
_site_lonlat['coj'] = (149.0708466, -31.2728196)
_site_lonlat['ogg'] = (-156.2589, 34.433161)

goodXTalkTargets = ['auto', 'GD 71', 'BD+284211', 'HZ 44', 'L745-46A', 'Feige 110', 'EGGR274']


def getAutoCandidate(context):
    if (context.site not in _site_lonlat):
        _logger.error("Site %s is not known. Giving up" % context.site)
        exit(1)

    site = ephem.Observer()
    lon, lat = _site_lonlat[context.site]
    site.lat = lat * math.pi / 180
    site.lon = lon * math.pi / 180
    site.date = ephem.Date(context.start + dt.timedelta(minutes=30))
    moon = ephem.Moon()
    moon.compute(site)

    print("Finding suitable star for site %s. Moon phase is  %i %%" % (context.site, moon.moon_phase * 100))

    for starcandidate in goodXTalkTargets:
        if 'auto' in starcandidate:
            continue
        radec = SkyCoord.from_name(starcandidate)
        s = ephem.FixedBody()
        s._ra = radec.ra.degree * math.pi / 180
        s._dec = radec.dec.degree * math.pi / 180
        s.compute(site)

        separation = (ephem.separation((moon.ra, moon.dec), (s.ra, s.dec)))
        alt = s.alt * 180 / math.pi
        separation = separation * 180 / math.pi

        altok = alt > 35
        sepok = separation > 30

        if (altok and sepok):
            print("\nViable star found: %s altitude % 4f moon separation % 4f" % (starcandidate, alt, separation))
            return starcandidate
        else:
            print("rejecting star %s - altitude ok: %s     moon separation ok: %s" % (starcandidate, altok, sepok))

    print("No viable star was found! full moon? giving up!")
    exit(1)


def createRequestsForStar(context):
    exposuretime = context.exptime
    overheadperexposure = 60 # readout plus fixed overhead
    telescopeslew = 120 + 90 + 60  # teelscope slew, acquire exposure + process, config change time
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

    pointing = {
        "type": "SP",
        "name": "Floyds test {}".format(context.targetname),
        "coord_type": "RD",
        "coord_sys": "ICRS",
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "pro_mot_ra": "0",
        "pro_mot_dec": "0",
        "parallax": "0.0000000",
        "ra": "%10f" % offsetPointing.ra.degree,
        "dec": "%10f" % offsetPointing.dec.degree,
    }

    print("Block for {} from {} to {}".format(pointing['name'], str(start), str(end)))

    # TODO: move to a better lcoation in code!
    agname = 'kb42' if 'ogg' in context.site else 'kb38'

    block = {"molecules": [],
             'start': start,
             'end': end,
             'site': context.site,
             'observatory': 'clma',
             'telescope': '2m0a',
             'instrument_class': '2M0-FLOYDS-SCICAM'.upper(),
             'priority': 31,
             "is_too": False,
             "max_airmass": "3.0000000",
             "min_lunar_dist": "30.0000000",
             }

    flat_molecule = {
        "spectra_slit": "slit_1.2as",
        "pointing": pointing,

        "tag_id": "LCOGT",
        "user_id": context.user,
        "prop_id": "LCOEngineering",
        "group": "Floyds test exposure",
        "exposure_count": 1,
        "bin_x": 1,
        "bin_y": 1,
        "inst_name": context.instrument,
        "priority": 1,
        "type": "LAMP_FLAT",
        "events": [],
        "ag_filter": "",
        "ag_exp_time": "10.0000000",
        "exposure_time": "40.0000000",
        "readout_mode": ""
    }

    spectrum_molecule = {
        "acquire_mode": "BRIGHTEST",
        "acquire_radius_arcsec": "5.00",
        "spectra_slit": "slit_1.2as",
        "pointing": pointing,
        "defocus": "0.0000000",
        "ag_name": agname,
        "ag_mode": "YES",
        "tag_id": "LCOGT",
        "user_id": context.user,
        "prop_id": "LCOEngineering",
        "group": "Floyds test exposure",
        "exposure_count": 1,
        "bin_x": 1,
        "bin_y": 1,
        "inst_name": context.instrument,
        "priority": 2,
        "type": "SPECTRUM",
        "events": [],
        "ag_filter": "",
        "ag_exp_time": "10.0000000",
        "exposure_time": exposuretime,
        "readout_mode": ""
    }

    arc_molecule = {
        "spectra_slit": "slit_1.2as",
        "pointing": pointing,
        "tag_id": "LCOGT",
        "user_id": context.user,
        "prop_id": "LCOEngineering-001",
        "group": "Floyds test exposure",
        "exposure_count": 1,
        "bin_x": 1,
        "bin_y": 1,
        "inst_name": context.instrument,
        "priority": 3,
        "type": "ARC",
        "events": [],
        "ag_filter": "",
        "ag_exp_time": "10.0000000",
        "exposure_time": "80.0000000",
        "readout_mode": ""
    }

    block['molecules'].append(flat_molecule)
    block['molecules'].append(spectrum_molecule)
    block['molecules'].append(arc_molecule)
    lastflat = copy.deepcopy (flat_molecule)
    lastflat['priority'] = 4
    block['molecules'].append(lastflat)

    _logger.debug(json.dumps(block, indent=4))
    if args.opt_confirmed:
        response = requests.post(LAKE_URL + '/blocks/', json=block)
        try:
            response.raise_for_status()
            _logger.info(
                'Submitted block with id: {0}. Check it at {1}/blocks/{0}'.format(response.json()['id'], LAKE_URL))
        except Exception:
            _logger.error(
                'Failed to submit block: error code {}: {}'.format(response.status_code, response.content))


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Submit an engineering Floyds observation.')

    parser.add_argument('--targetname', default='auto', type=str,
                        help='Name of star for Floyds test observation. if none is given, or auto, Software will try to'
                             ' find a cataloged flux stadnard star. Name must resolve via Simbad; if resolve failes,'
                             ' program will exit.')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--site', choices=['ogg', 'coj'], required=True,  help="To which site to submit")


    parser.add_argument('--exp-cnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, default=150)

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
            _logger.error("Invalidt start time argument: ", args.start)
            exit(1)

    if ('auto' in args.targetname):
        # automatically find the best target
        args.targetname = getAutoCandidate(args)
        pass

    try:
        _logger.debug("Resolving target name")
        args.radec = SkyCoord.from_name(args.targetname)
    except:
        print("Resolving target name failed, giving up")
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.targetname, args.radec.ra, args.radec.dec))
    return args


if __name__ == '__main__':
    args = parseCommandLine()

    createRequestsForStar(args)
