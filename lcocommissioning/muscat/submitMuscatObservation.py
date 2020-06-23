import argparse
import datetime as dt
import json
import logging

from astropy.coordinates import SkyCoord

import lcocommissioning.common.common as common

_log = logging.getLogger(__name__)


def createMuscatRequest(args):
    configurations = {
        'configurations': [{
            'type': None,
            'instrument_type': '2M0-SCICAM-MUSCAT',
            'guiding_config': {},
            'instrument_configs': [{
                'exposure_count': 1,
                'optical_elements': {
                    'diffuser_g_position': 'in',
                    'diffuser_r_position': 'in',
                    'diffuser_i_position': 'in',
                    'diffuser_z_position': 'in',
                },
                'extra_params': {
                    'exposure_time_g': args.exp_times[0],
                    'exposure_time_r': args.exp_times[1],
                    'exposure_time_i': args.exp_times[2],
                    'exposure_time_z': args.exp_times[3],
                    'exposure_mode': args.exp_mode,
                }
            }, ]
        }, ]
    }

    if args.exp_cnt:
        configurations['configurations'][0]['type'] = 'EXPOSE'
        configurations['configurations'][0]['instrument_configs'][0]['exposure_count'] = args.exp_cnt

    elif args.filltime:
        configurations['configurations'][0]['type'] = 'REPEAT_EXPOSE'
        configurations['configurations'][0]['repeat_duration'] = args.filltime

    return configurations


def createRequest(args):
    pass


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='MuSCAT @ LCO engineering commissioning submisison tool')

    parser.add_argument('--targetname', default='auto', type=str,
                        help='Name of object to observe; will beresolved via simbad. Can be coordinates in the form of Jhhmmss+ddmmss')
    parser.add_argument('--title', default="Muscat commissioning", help="Descriptive title for observation request")

    parser.add_argument('--site', default='ogg', choices=common.lco_2meter_sites,
                        help="To which site to submit")

    # parser.add_argument('--readmode', choices=common.lco_muscat_readmodes, default=common.lco_muscat_readmodes[0])

    parser.add_argument('--start', default=None,
                        help="Time to start observation. If not given, defaults to \"NOW\". Specify as YYYYMMDD HH:MM")

    parser.add_argument('--schedule-window', default=2, type=float, help="How long after start should request be schedulable?")

    # parser.add_argument('--dither', action='store_true', help='Dither the exposure in a 5 point pattern.')

    parser.add_argument('--defocus', type=float, default=0.0, help="Amount to defocus star.")

    parser.add_argument('--exp-times', nargs=4, type=float, default=[2, 4, 6, 12],
                        help='List of exposure times in g r i z, e.g., "1.4 1.6 2.0 5')

    parser.add_argument('--ipp', type=float, default=1.0, help="ipp value")


    parser.add_argument('--offsetRA', default=0, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=0, help="Extra pointing offset to apply Dec")
    parser.add_argument('--pp', default=[0., 0.], nargs=2, type=float, help="Proper motion, mas/yr")


    parser.add_argument('--exp-mode', default='SYNCHRONOUS', choices=['SYNCHRONOUS', 'ASYNCHRONOUS'], required=False)

    repeatgroup = parser.add_mutually_exclusive_group()
    repeatgroup.add_argument('--exp-cnt', type=int, help="How often to reapeat each exposure")
    repeatgroup.add_argument('--filltime', type=float, help="How long to repeat Muscat exposures (seconds)")


    parser.add_argument('--scheduler', action='store_true',
                        help='If set, submit to scheduler instead of tryuing a direct submission.')
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
            exit(1)

    try:
        _log.debug("Resolving target name")
        args.radec = SkyCoord.from_name(args.targetname, parse=True)
    except:
        print("Resolving target name failed, giving up")
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.targetname, args.radec.ra, args.radec.dec))

    if not (args.exp_cnt  or args.filltime):
        print ("No exposure mode choosen, defaulting to EXPOSE")
        args.exp_cnt = 1

    return args


def main():
    args = parseCommandLine()

    muscat = createMuscatRequest(args)
    print(json.dumps(muscat, indent=1))  # if args.scheduler:
    #     create_request_for_star_scheduler(args)
    # else:
    #     createRequestsForStar_pond(args)


if __name__ == '__main__':
    main()
