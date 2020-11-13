import argparse
import json
import logging

import astropy
from astropy.coordinates import SkyCoord
import datetime as dt
import lcocommissioning.common.common as common

_logger = logging.getLogger(__name__)

defaultconstraints = {"max_airmass": 2.5,
                      "min_lunar_distance": 45.0, }


def createNRESRequestsConfiguration(context):
    configuration = {
        'type': 'NRES_SPECTRUM',
        'instrument_type': '1M0-NRES-SCICAM',
        'guiding_config': {'mode': 'ON'},
        'acquisition_config': {'mode': 'WCS' if context.forcewcs else 'BRIGHTEST'},
        'instrument_configs': [{
            'exposure_time': context.exptime,
            'exposure_count': context.expcnt,
        },]
    }
    return configuration


def createRequest(args):
    requestgroup = {"name": f'NRES engineering {args.targetname}',
                    "proposal": "ENG2017AB-001",
                    "ipp_value": args.ipp,
                    "operator": "SINGLE",  # "MANY" if args.dither else "SINGLE",
                    "observation_type": "NORMAL",
                    "requests": []
                    }

    absolutestart = args.start
    windowend = args.start + dt.timedelta(hours=args.schedule_window)

    location = {
        'telescope_class': '1m0',
    }
    if args.site is not None:
        location['site'] = args.site

    request = {'configurations': [],
               'windows': [{"start": absolutestart.isoformat(), "end": windowend.isoformat()}, ],
               'location': location}

    nresconfiguration = createNRESRequestsConfiguration(args)

    pm_ra = args.pm[0]
    pm_dec = args.pm[1]
    target = {
        "type": "ICRS",
        "name": f"NRES Commissioning {args.targetname}",
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "ra": "%10f" % args.radec.ra.degree,
        "dec": "%10f" % args.radec.dec.degree,
        "proper_motion_ra": pm_ra,
        "proper_motion_dec": pm_dec,
    }

    nresconfiguration['target'] = target
    nresconfiguration['constraints'] = common.default_constraints
    request['configurations'].append(nresconfiguration)
    requestgroup['requests'].append(request)
    return requestgroup


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
    parser.add_argument('--pm', type=float, nargs=2, help="proper motion RA DEC in marcsec / year")

    parser.add_argument('--start', default=None,
                        help="When to start NRES observation. If not given, defaults to \"NOW\"")
    parser.add_argument('--schedule-window', default=3, type=float,
                        help="How long after start should request be schedulable?")

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

    if args.targetname is None:
        args.targetname = common.get_auto_target(common.goodXTalkTargets, site=args.site, starttime=args.start)
        _logger.info(f"Auto selecting target {args.targetname}")

    if args.targetname is None:
        _logger.error("No target given, giving up.")
        exit(1)

    astropy.coordinates.name_resolve.sesame_database.set("simbad")
    try:
        _logger.debug("Resolving target name")

        args.radec = SkyCoord.from_name(args.targetname, parse=True)
    except Exception as e:
        print("Resolving target name failed, giving up {}".format(e))
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.targetname, args.radec.ra, args.radec.dec))
    return args


def main():
    args = parseCommandLine()
    requstgroup = createRequest(args)
    common.send_request_to_portal(requstgroup, args.opt_confirmed)
    exit(0)


if __name__ == '__main__':
    main()
