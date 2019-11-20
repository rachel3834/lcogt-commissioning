import argparse
import json
import logging
import ephem
import math
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import datetime as dt
import lcocommissioning.common.common as common

_log = logging.getLogger(__name__)

sinistro_1m_quadrant_offsets = {0: [-60, 60],
                                1: [60, 60],
                                2: [60, -60],
                                3: [-60, -60]}

no_dither = {0: [0, 0], }


def getRADecForQuadrant(starcoo, quadrant, extraoffsetra=0, extraoffsetDec=0):
    dra = Angle(sinistro_1m_quadrant_offsets[quadrant][0], unit=u.arcsec)
    ddec = Angle(sinistro_1m_quadrant_offsets[quadrant][1], unit=u.arcsec)

    return SkyCoord(starcoo.ra - dra + Angle(extraoffsetra, unit=u.arcsec),
                    starcoo.dec - ddec + Angle(extraoffsetDec, unit=u.arcsec))


def create_request_for_star_scheduler(context):
    absolutestart = context.start
    windowend = context.start + dt.timedelta(hours=context.schedule_window)

    location = {'telescope': '1m0a',
                'telescope_class': '1m0',
                'site': context.site, }

    if context.dome != "None":
        location['enclosure'] = context.dome

    if context.instrument != "None":
        location['instrument'] = context.instrument

    requestgroup = {"name": context.title,
                    "proposal": "LCOEngineering",
                    "ipp_value": context.ipp,
                    "operator": "SINGLE" if context.nodither else "MANY",
                    "observation_type": "NORMAL",
                    "requests": []
                    }

    offsets = sinistro_1m_quadrant_offsets
    if context.nodither:
        offsets = no_dither

    for quadrant in offsets:
        request = {'configurations': None,
                    'windows': [{"start": str(absolutestart), "end": str(windowend)}, ],
                    'location': location}

        offsetPointing = getRADecForQuadrant(context.radec, quadrant, context.offsetRA, context.offsetDec)
        target = {
            "type": "ICRS",
            "name": "{} {}".format(context.title, context.targetname),
            "epoch": "2000.0000000",
            "equinox": "2000.0000000",
            "ra": "%10f" % offsetPointing.ra.degree,
            "dec": "%10f" % offsetPointing.dec.degree,
        }

        configurations = []
        for exptime in context.exp_times:
            configuration = {
            'type': 'EXPOSE',
            'instrument_type': '1M0-SCICAM-SINISTRO',
            'target': target,
            'constraints': common.default_constraints,
            'acquisition_config': {},
            'guiding_config': {},
            'instrument_configs': [],
            }

            configuration['instrument_configs'].append(
                {
                    'exposure_time': exptime,
                    'exposure_count': context.exp_cnt,
                    'mode': context.readmode,
                    'optical_elements': {
                        'filter': context.filter
                    },
                    "extra_params": {
                        "defocus": min(5, context.defocus)  # scheduler doe snot allow defocussing more than 3mm FP.
                    }
                })

            configurations.append (configuration)
        request['configurations']=configurations
        requestgroup['requests'].append (request)

    common.send_request_to_portal(requestgroup, context.opt_confirmed)


def createRequestsForStar_pond(context):
    timePerQuadrant = 8  # in minutes
    absolutestart = context.start

    # create one block per quadrant
    offsets = sinistro_1m_quadrant_offsets
    if context.nodither:
        offsets = no_dither
    for quadrant in offsets:

        start = absolutestart + dt.timedelta(minutes=(quadrant) * timePerQuadrant)
        end = start + dt.timedelta(minutes=timePerQuadrant)
        start = str(start).replace(' ', 'T')
        end = str(end).replace(' ', 'T')

        print("Block for %s Q %d from %s to %s" % (context.targetname, quadrant, str(start), str(end)))

        block_params = {
            "molecules": [],
            'start': start,
            'end': end,
            'site': context.site,
            'observatory': context.dome,
            'telescope': context.telescope,
            'instrument_class': '1m0-SciCam-Sinistro'.upper(),
            'priority': 30,
        }

        offsetPointing = getRADecForQuadrant(context.radec, quadrant, context.offsetRA, context.offsetDec)

        for exptime in context.exp_times:
            moleculeargs = {
                'inst_name': context.instrument,
                'bin': 2,
                'exposure_time': exptime,
                "readout_mode": context.readmode,
                'exposure_count': context.exp_cnt,
                'bin_x': 2,
                'bin_y': 2,
                #         'ag_mode': 'ON',

                'filter': context.filter,
                'pointing': {"type": "SP",
                             "name": "%s x talk q %d" % (context.targetname, quadrant),
                             "coord_type": "RD",
                             "coord_sys": "ICRS",
                             "epoch": "2000",
                             "equinox": "2000",
                             "ra": "%10f" % offsetPointing.ra.degree,
                             "dec": "%7f" % offsetPointing.dec.degree,
                             "proper_motion_ra": "%7.3f" % context.pp[0],
                             "proper_motion_dec": "%7.3f" % context.pp[1],
                             },

                'group': 'Sinistro x talk commissioning',
                'user_id': context.user,
                'prop_id': 'LCOEngineering',
                'defocus': context.defocus,
                'type': 'EXPOSE'
            }

            block_params['molecules'].append(moleculeargs)

        common.send_to_lake(block_params, context.opt_confirmed)


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='X-Talk calibration submission tool\nSubmit to LAKE the request to observe a bright star, '
                    'defocussed, at 1,3,6,12 sec exposure time, on each quadrant. Useful when commissioing a camera '
                    'that is not available via scheduler yet.')

    parser.add_argument('--site', required=True, choices=common.lco_1meter_sites,
                        help="To which site to submit")
    parser.add_argument('--dome', required=True, choices=['doma', 'domb', 'domc', 'None'],
                        help="To which enclosure to submit")
    parser.add_argument('--telescope', default='1m0a')
    parser.add_argument('--instrument', required=True,
                        choices=common.lco_sinistro1m_cameras.extend('None'),
                        help="To which instrument to submit")
    parser.add_argument('--readmode', choices=common.archon_readout_modes, default=common.archon_readout_modes[0])
    parser.add_argument('--targetname', default='auto', type=str,
                        help='Name of star for X talk measurement. Will be resolved via simbad. If resolve failes, '
                             'program will exit.\n If name is auto, which is te default, a viable target will be choosen for you.')
    parser.add_argument('--title', default="Sinsitro commissioning: X talk")
    parser.add_argument('--start', default=None,
                        help="Time to start x-talk calibration. If not given, defaults to \"NOW\". Specify as YYYYMMDD HH:MM")

    parser.add_argument('--nodither', action='store_true',
                        help='Do not dither the exposure')
    parser.add_argument('--defocus', type=float, default=6.0, help="Amount to defocus star.")
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")
    parser.add_argument('--exp-times', nargs="*", type=float, default=[2, 4, 6, 12], help="List of exposure times")
    parser.add_argument('--exp-cnt', type=int, default=1, help="How often to reapeat each exposure")
    parser.add_argument('--ipp', type=float, default=1.0, help="ipp value")
    parser.add_argument('--filter', type=str, default='rp', help="Filter")
    parser.add_argument('--offsetRA', default=0, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=0, help="Extra pointing offset to apply Dec")
    parser.add_argument('--pp', default=[0., 0.], nargs=2, type=float, help="Proper motion, mas/yr")
    parser.add_argument('--schedule-window', default=2, type=float)
    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, block will be submitted. If omitted, nothing will be submitted.')
    parser.add_argument('--scheduler', action='store_true',
                        help='If set, submit to scheduler instead of Lake.')

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


if __name__ == '__main__':
    main()
