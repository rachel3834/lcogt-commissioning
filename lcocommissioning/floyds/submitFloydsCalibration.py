import argparse
import logging
import datetime as dt

import lcocommissioning.common.common as common

_logger = logging.getLogger(__name__)


def makeDirectFlatFieldsubmission(args):
    start = str(args.start+ dt.timedelta(seconds=30)).replace(' ', 'T')
    end = args.start + dt.timedelta(seconds= args.expcnt * (args.exptime + 30))
    end = str(end).replace(' ','T')
    observation = {
        'name': 'Floyds direct flat field testing',
        'proposal': 'LCOEngineering',
        'site': args.site,
        'enclosure': 'clma',
        'telescope': '2m0a' ,
        'start': start,
        'end': end,
        'priority': 60,
        'request': {
            'configurations': [{
                'type': 'LAMP_FLAT',
                'instrument_type': '2M0-FLOYDS-SCICAM'.upper(),
                'instrument_name': args.instrument,
                'instrument_configs': [{
                    'exposure_count': args.expcnt,
                    'exposure_time': args.exptime,
                    'bin_x': 1,
                    'bin_y': 1,
                    'optical_elements': {
                        'slit': args.slit
                    }

                 }],

            'target': {
                'type': 'HOUR_ANGLE',
                'name': 'Fake target',
                'hour_angle': 0,
                'dec': 0
            },

            'acquisition_config': {},
            'guiding_config': {},
            'constraints': {}
            }]
        }
    }
    _logger.info (f"Start time is {start}")
    common.submit_observation(observation, args.opt_confirmed)



def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='Submit an engineering Floyds Flat (Tungsten lamp) observation.')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('--site', choices=common.lco_2meter_sites, required=True,  help="To which site to submit")
    parser.add_argument('--exp-cnt', type=int, dest="expcnt", default=1)
    parser.add_argument('--exptime', type=float, default=150)
    parser.add_argument('--slit', type=str, default="slit_1.2as", choices=['slit_1.2as','slit_2.0as','slit_6.0as'])

    parser.add_argument('--start', default=None,
                        help="When to start Floyds observation. If not given, defaults to \"NOW\"")
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")

    # Per default, do not be on chip gap!

    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, block will be submitted.')

    parser.add_argument('--loglevel', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
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

    return args


def main():
    args = parseCommandLine()

    makeDirectFlatFieldsubmission(args)


if __name__ == '__main__':
    main()
