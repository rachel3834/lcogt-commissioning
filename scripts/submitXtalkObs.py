import argparse
import logging
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import datetime as dt
from lcogtpond import block, molecule, pointing


_logger = logging.getLogger(__name__)

quadrantOffsets = {0: [-450, 450],
                   1: [450, 450],
                   2: [450, -450],
                   3: [-450, -450]}


def getRADecForQuadrant(starcoo, quadrant):
    dra = Angle(quadrantOffsets[quadrant][0], unit=u.arcsec)
    ddec = Angle(quadrantOffsets[quadrant][1], unit=u.arcsec)

    return SkyCoord(starcoo.ra + dra, starcoo.dec + ddec)




def createRequestsForStar(context):

    if context.start is None:
        start = dt.datetime.utcnow()
    else:
        try:
            start = dt.datetime.strptime(context.start, "%Y%m%d %H:%M")
        except ValueError:
            _logger.error ("Invalidt start time argument: ", context.start)
            exit(1)

    end = start + dt.timedelta(hours=1)

    print ("Block for %s from %s to %s" % (context.name, start, end))

    block_params = {
        'start'       : start,
        'end'         : end,
        'site'        : context.site,
        'observatory' : context.dome,
        'telescope'   : context.telescope,
        'instrument_class': '1m0-SciCam-Sinistro'.upper(),
        'priority'    : 30,
    }
    my_block = block.Block.build(**block_params)




    for quadrant in quadrantOffsets:

        offsetPointing = getRADecForQuadrant(context.radec, quadrant)

        for exptime in [1,3,6,12]:

            moleculeargs = {
                'inst_name' : context.instrument,
                'bin' : 1,
                'exp_time': exptime,
                'exp_cnt': 1,
                'filters' : 'rp',
                'pointing': pointing.sidereal(name="%s x talk q %d" % (context.name, quadrant),
                                              coord=pointing.ra_dec(ra=offsetPointing.ra.degree, dec=offsetPointing.dec.degree)),

                'group': 'Sinistro x talk commissioing',
                'user': context.user,
                'proposal' : 'calibration',
                'ag_mode': 0,

                'defocus': context.defocus

            }

            my_molecule = molecule.Expose.build(**moleculeargs)

            my_block.add_molecule (my_molecule)



    if args.opt_confirmed:
        print ("Save block ...")
        my_block.save()



def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='X-Talk calibration submission tool\nSubmit to POND the request to observe a bright star, defocussed, at 1,3,6,12 sec exposure time, on each quadrant.')

    parser.add_argument('--name', required=True,  type=str, help='Name of star for X talk measurement. Will be resolved via simbad. If resolve failes, program will exit.\n future version will automatically select a star based on time of observation.  ')
    parser.add_argument('--defocus', type=float, default=6.0, help="Amount to defocus star.")
    parser.add_argument('--site', default='coj', choices=['lsc', 'cpt', 'coj', 'elp'], help="To which site to submit")
    parser.add_argument('--dome', default='doma',choices=['doma', 'domb', 'domc'], help="To which enclosure to submit")
    parser.add_argument('--telescope', default='1m0a')
    parser.add_argument('--instrument', default='fl12', choices=['fl03', 'fl04', 'fl05', 'fl12', 'fl15', 'fl16', ] ,help="To which instrumetn to submit")
    parser.add_argument('--start', default=None, help="When to start x-talk calibration. If not given, defaults to \"NOW\"")
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")
    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, block will be submitted.')

    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')

    try:
        _logger.debug("Resolving target name")
        args.radec = SkyCoord.from_name(args.name)
    except:
        print("Resolving target name failed, giving up")
        exit(1)

    print ("Resolved target %s at corodinates %s %s" % (args.name, args.radec.ra, args.radec.dec))
    return args


if __name__ == '__main__':
    args = parseCommandLine()

    createRequestsForStar(args)



