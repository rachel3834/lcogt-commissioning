import argparse
import datetime as dt
import json
import logging

import numpy as np
from astropy.coordinates import SkyCoord

import lcocommissioning.common.common as common

_log = logging.getLogger(__name__)


def createMuscatRequestConfiguration(args):
    configuration = {
            'type': None,
            'instrument_type': '2M0-SCICAM-MUSCAT',
            'guiding_config': {},
            'acquisition_config': {},
            'instrument_configs': [{
                'exposure_count': 1,
                'optical_elements': {
                    'diffuser_g_position': 'out',
                    'diffuser_r_position': 'out',
                    'diffuser_i_position': 'out',
                    'diffuser_z_position': 'out',
                },
                'extra_params': {
                    'exposure_time_g': args.exp_times[0],
                    'exposure_time_r': args.exp_times[1],
                    'exposure_time_i': args.exp_times[2],
                    'exposure_time_z': args.exp_times[3],
                    'exposure_mode': args.exp_mode,
                }
            }, ]
        }

    if args.exp_cnt:
        configuration['type'] = 'EXPOSE'
        configuration['instrument_configs'][0]['exposure_count'] = args.exp_cnt

    elif args.filltime:
        configuration['type'] = 'REPEAT_EXPOSE'
        configuration['repeat_duration'] = args.filltime

    return configuration


def createRequest(args):
    requestgroup = {"name": args.title,
                    "proposal": "MuSCAT Commissioning",
                    "ipp_value": args.ipp,
                    "operator": "SINGLE", # "MANY" if args.dither else "SINGLE",
                    "observation_type": "NORMAL",
                    "requests": []
                    }

    absolutestart = args.start
    windowend = args.start + dt.timedelta(hours=args.schedule_window)

    location = {'telescope': '2m0a',
                'telescope_class': '2m0',
                'enclosure' : 'clma',
                'site': args.site, }

    request = {'configurations': [],
               'windows': [{"start": absolutestart.isoformat(), "end": windowend.isoformat()}, ],
               'location': location}

    muscatconfiguration = createMuscatRequestConfiguration(args)

    target = {
        "type": "ICRS",
        "name": "{} {}".format(args.title, args.targetname),
        "epoch": "2000.0000000",
        "equinox": "2000.0000000",
        "ra": "%10f" % args.radec.ra.degree,
        "dec": "%10f" % args.radec.dec.degree,
    }
    muscatconfiguration['target'] = target
    muscatconfiguration['constraints'] = common.default_constraints
    request['configurations'].append (muscatconfiguration)

    requestgroup['requests'].append (request)

    return requestgroup


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='MuSCAT @ LCO engineering commissioning submisison tool')

    parser.add_argument('--targetname', default='auto', type=str,
                        help='Name of object to observe; will beresolved via simbad. Can be coordinates in the form of Jhhmmss+ddmmss')
    parser.add_argument('--title', default="Muscat commissioning", help="Descriptive title for observation request")

    parser.add_argument('--site', default='ogg', choices=common.lco_2meter_sites,
                        help="To which site to submit")

    parser.add_argument('--start', default=None,
                        help="Time to start observation. If not given, defaults to \"NOW\". Specify as YYYYMMDD HH:MM")

    parser.add_argument('--schedule-window', default=2, type=float,
                        help="How long after start should request be schedulable?")

    # parser.add_argument('--dither', action='store_true', help='Dither the exposure in a 5 point pattern.')

    parser.add_argument('--defocus', type=float, default=0.0, help="Amount to defocus star.")

    parser.add_argument('--exp-times', nargs=4, type=float, default=[10, 10, 10, 10],
                        help='List of exposure times in g r i z, e.g., "1.4 1.6 2.0 5')

    parser.add_argument('--ipp', type=float, default=1.0, help="ipp value")

    parser.add_argument('--offsetRA', default=0, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=0, help="Extra pointing offset to apply Dec")
    parser.add_argument('--pp', default=[0., 0.], nargs=2, type=float, help="Proper motion, mas/yr")

    parser.add_argument('--exp-mode', default='SYNCHRONOUS', choices=['SYNCHRONOUS', 'ASYNCHRONOUS'], required=False)

    repeatgroup = parser.add_mutually_exclusive_group()
    repeatgroup.add_argument('--exp-cnt', type=int, help="How often to repeat each exposure")
    repeatgroup.add_argument('--filltime', type=float, help="How long to repeat Muscat exposures (seconds)")

    parser.add_argument('--scheduler', action='store_true',
                        help='If set, submit to scheduler instead of trying a direct submission.')
    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, observation will be submitted. If omitted, nothing will be submitted.')

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
            _log.error("Invalid start time argument: ", args.start)
            exit(1)

    if ('auto' in args.targetname):
        # automatically find the best target
        args.targetname = common.get_auto_target(common.goodXTalkTargets, args.site, args.start, moonseparation=40)
        if args.targetname is None:
            _log.error ("Could not find a suitable auto target. Exiting.")
            exit(1)

    try:
        _log.debug("Resolving target name")
        args.radec = SkyCoord.from_name(args.targetname, parse=True)
    except:
        print("Resolving target name failed, giving up")
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.targetname, args.radec.ra, args.radec.dec))

    if not (args.exp_cnt or args.filltime):
        print("No exposure mode choosen, defaulting to EXPOSE")
        args.exp_cnt = 1

    return args


def ammend_request_for_direct_submission(muscat, args):
    """ Need to fill in some extra fields for a direct submission"""

    # Calculate the4 end tim efor this reuquest
    if args.filltime:
        endtime = args.start + dt.timedelta(seconds=args.filltime)
    if args.exp_cnt:
        endtime = args.start + dt.timedelta (seconds = args.exp_cnt * (float (np.max (args.exp_times) + 10)))

    data = {
        'name': args.title,
        'proposal': 'LCOEngineering',
        'site': args.site,
        'enclosure': 'clma',
        'telescope': '2m0a',
        'start': args.start.isoformat(),
        'end': endtime.isoformat(),
        'request': {
            'acceptability_threshold' : 90,
            'configurations': muscat['requests'][0]['configurations']
        }

    }
    return data


def main():
    args = parseCommandLine()

    muscat = createRequest(args)

    if args.scheduler:
        _log.info ("Submitting to scheduler")
        _log.info(json.dumps(muscat, indent=2))
        common.send_request_to_portal(muscat, args.opt_confirmed)
    else:
        _log.info ("Attempting direct submission")
        muscat = ammend_request_for_direct_submission (muscat, args)
        _log.info(json.dumps(muscat, indent=2))
        common.submit_observation(muscat, args.opt_confirmed)


# else:
    #     createRequestsForStar_pond(args)


if __name__ == '__main__':
    main()
