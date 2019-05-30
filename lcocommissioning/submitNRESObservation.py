import argparse
import json
import logging

import astropy
from astropy.coordinates import SkyCoord
import datetime as dt
import lcocommissioning.common.common as common

_logger = logging.getLogger(__name__)

defaultconstraints = {"max_airmass": 2.0,
                      "min_lunar_distance": 30.0, }


def createRequestsForStar(context):
    start = context.start
    end = start + dt.timedelta(minutes=60 * context.window)
    start = str(start).replace(' ', 'T')
    end = str(end).replace(' ', 'T')
    targetPointing = SkyCoord(context.radec.ra, context.radec.dec)

    pointing = {
        "type": "SIDEREAL",
        "name": "NRES commissioning {}".format(context.targetname),
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "pro_mot_ra": "0",
        "pro_mot_dec": "0",
        "parallax": "0.0000000",
        "ra": "%10f" % targetPointing.ra.degree,
        "dec": "%10f" % targetPointing.dec.degree,
    }

    nres_molecule = {
        "type": "NRES_SPECTRUM",
        "ag_mode": "ON",
        "ag_strategy": "",
        "instrument_name": "1M0-NRES-SCICAM",
        "acquire_mode": "BRIGHTEST",
        "acquire_radius_arcsec": 5.0,
        "acquire_exp_time": None,
        "expmeter_mode": "OFF",
        "exposure_time": context.exptime,
        "exposure_count": context.expcnt,
        "bin_x": 1,
        "bin_y": 1,
    }

    if context.forcewcs:
        nres_molecule['ag_strategy'] = 'wcs'
      #  nres_molecule['ag_strategy'] = None # test out None acquistion strategy as per Steve's recommendation.
        nres_molecule['acquire_strategy'] = 'astrometry'

    userrequest = {"group_id": "NRES test observation",
                   "proposal": context.proposalid,
                   "ipp_value": 1.1,
                   "operator": "SINGLE",
                   "observation_type": "NORMAL",
                   "requests": [
                       {"acceptability_threshold": 90,
                        "target": pointing,
                        "molecules": [nres_molecule, ],
                        "windows": [{"start": str(start), "end": str(end)}, ],
                        "location": {'telescope_class': '1m0'},
                        "constraints": defaultconstraints,
                        },
                   ]}

    if context.site is not None:
        userrequest['requests'][0]['location']['site'] = context.site
    _logger.debug(json.dumps(userrequest, indent=4))

    common.send_to_scheduler(userrequest, context.opt_confirmed)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Submit an engineering NRES observation to SCHEDULER.')

    parser.add_argument("--proposalid", default="ENG2017AB-001")

    parser.add_argument('--targetname', type=str,
                        help='Name of star for NRES test observation. if none is given, or auto, Software will try to'
                             ' find a cataloged flux stadnard star. Name must resolve via Simbad; if resolve failes,'
                             ' program will exit.')

    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument('--site', choices=['lsc', 'elp', 'cpt', 'tlv'], default=None,
                               help="To which site to submit")
    parser.add_argument('--exp-cnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, default=120)
    parser.add_argument('--forcewcs', action='store_true',
                        help='Force WCSW based acquistion')
    parser.add_argument('--window', default=3, type=int, help="scheduling window length")
    parser.add_argument('--ipp', default=1.0, help="IPP priority for block")

    parser.add_argument('--start', default=None,
                        help="When to start Floyds observation. If not given, defaults to \"NOW\"")
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")
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
            _logger.error("Invalid start time argument: ", args.start)
            exit(1)
    astropy.coordinates.name_resolve.sesame_database.set("simbad")
    try:
        _logger.debug("Resolving target name")

        args.radec = SkyCoord.from_name(args.targetname)
    except Exception as e:
        print("Resolving target name failed, giving up {}".format (e))
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.targetname, args.radec.ra, args.radec.dec))
    return args


def main():
    args = parseCommandLine()
    createRequestsForStar(args)

if __name__ == '__main__':
    main()
