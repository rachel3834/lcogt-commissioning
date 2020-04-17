import argparse
import lcocommissioning.common.common as common


def parseCommandLine():

    parser = argparse.ArgumentParser(
        description='MuSCAT @ LCO engineering commissioning submisison tool')

    parser.add_argument('--site', default='ogg', choices=common.lco_2meter_sites,
                        help="To which site to submit")
    parser.add_argument('--enclosure', default='clma', choices=['clma'],
                        help="To which enclosure to submit")
    parser.add_argument('--telescope', default='2m0a')
    parser.add_argument('--instrument', required=True,
                        choices=common.lco_muscat_instruments.extend('None'),
                        help="To which instrument to submit")
    parser.add_argument('--readmode', choices=common.lco_muscat_readmodes, default=common.lco_muscat_readmodes[0])
    parser.add_argument('--targetname', default='auto', type=str,
                        help='Name of object to observe; will beresolved via simbad. Can be coordinates in the form of Jhhmmss+ddmmss')
    parser.add_argument('--title', default="Muscat commissioning")
    parser.add_argument('--start', default=None,
                        help="Time to start observation. If not given, defaults to \"NOW\". Specify as YYYYMMDD HH:MM")

    parser.add_argument('--dither', action='store_true',
                        help='Dither the exposure')
    parser.add_argument('--defocus', type=float, default=0.0, help="Amount to defocus star.")

    parser.add_argument('--exp-times', nargs="*", type=float, default=[2, 4, 6, 12], help="List of exposure times")
    parser.add_argument('--exp-cnt', type=int, default=1, help="How often to reapeat each exposure")
    parser.add_argument('--ipp', type=float, default=1.0, help="ipp value")
    parser.add_argument('--filter', type=str, default='rp', help="Filter")
    parser.add_argument('--offsetRA', default=0, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=0, help="Extra pointing offset to apply Dec")
    parser.add_argument('--pp', default=[0., 0.], nargs=2, type=float, help="Proper motion, mas/yr")
    parser.add_argument('--schedule-window', default=2, type=float)
    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, observation will be submitted. If omitted, nothing will be submitted.')
    parser.add_argument('--scheduler', action='store_true',
                        help='If set, submit to scheduler instead of directly.')

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
    return args




def main():
    args = parseCommandLine()
    if args.scheduler:
        create_request_for_star_scheduler(args)
    else:
        createRequestsForStar_pond(args)
