import argparse
import logging
import ephem
import math
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import datetime as dt
from lcogtpond import block, molecule, pointing

_logger = logging.getLogger(__name__)

quadrantOffsets = {0: [-450, 450],
                   1: [450, 450],
                   2: [450, -450],
                   3: [-450, -450]}

_site_lonlat = {}
_site_lonlat['bpl'] = (-119.863103, 34.433161)
_site_lonlat['coj'] = (149.0708466, -31.2728196)
_site_lonlat['cpt'] = (20.8124, -32.3826)
_site_lonlat['elp'] = (-104.015173, 30.679833)
_site_lonlat['lsc'] = (-70.8049, -30.1673666667)
_site_lonlat['ogg'] = (-156.2589, 34.433161)
_site_lonlat['sqa'] = (-120.04222167, 34.691453333)
_site_lonlat['tfn'] = (-16.511544, 28.300433)

goodXTalkTargets = ['auto', '91 Aqr', 'HD30562', '15 Sex']


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

        altok = alt > 40
        sepok = separation > 30

        if (altok and sepok):
            print("\nViable star found: %s altitude % 4f moon separation % 4f" % (starcandidate, alt, separation))
            return starcandidate
        else:
            print("rejecting star %s - altitude ok: %s     moon separation ok: %s" % (starcandidate, altok, sepok))

    print("No viable star was found! full moon? giving up!")
    exit(1)


def getRADecForQuadrant(starcoo, quadrant, extraoffsetra=0, extraoffsetDec=0):
    dra = Angle(quadrantOffsets[quadrant][0], unit=u.arcsec)
    ddec = Angle(quadrantOffsets[quadrant][1], unit=u.arcsec)

    return SkyCoord(starcoo.ra - dra + Angle(extraoffsetra, unit=u.arcsec),
                    starcoo.dec - ddec + Angle(extraoffsetDec, unit=u.arcsec))


def createRequestsForStar(context):
    timePerQuadrant = 8  # in minutes
    absolutestart = context.start

    # create one block per quadrant
    for quadrant in quadrantOffsets:

        start = absolutestart + dt.timedelta(minutes=(quadrant) * timePerQuadrant)
        end = start + dt.timedelta(minutes=timePerQuadrant)

        print("Block for %s Q %d from %s to %s" % (context.name, quadrant, start, end))

        block_params = {
            'start': start,
            'end': end,
            'site': context.site,
            'observatory': context.dome,
            'telescope': context.telescope,
            'instrument_class': '1m0-SciCam-Sinistro'.upper(),
            'priority': 40,
        }

        my_block = block.Block.build(**block_params)

        offsetPointing = getRADecForQuadrant(context.radec, quadrant, context.offsetRA, context.offsetDec)

        for exptime in [2, 4, 6, 12]:
            moleculeargs = {
                'inst_name': context.instrument,
                'bin': 1,
                'exp_time': exptime,
                'exp_cnt': 1,
                'filters': 'rp',
                'pointing': pointing.sidereal(name="%s x talk q %d" % (context.name, quadrant),
                                              coord=pointing.ra_dec(ra=offsetPointing.ra.degree,
                                                                    dec=offsetPointing.dec.degree)),

                'group': 'Sinistro x talk commissioning',
                'user': context.user,
                'proposal': 'calibration',
                'defocus': context.defocus
            }

            my_molecule = molecule.Expose.build(**moleculeargs)
            my_block.add_molecule(my_molecule)

        if args.opt_confirmed:
            print("Save block ...")
            my_block.save()


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='X-Talk calibration submission tool\nSubmit to POND the request to observe a bright star, defocussed, at 1,3,6,12 sec exposure time, on each quadrant.')

    parser.add_argument('--name', default='auto', type=str, choices=goodXTalkTargets,
                        help='Name of star for X talk measurement. Will be resolved via simbad. If resolve failes, program will exit.\n future version will automatically select a star based on time of observation.  ')
    parser.add_argument('--defocus', type=float, default=6.0, help="Amount to defocus star.")
    parser.add_argument('--site', default='coj', choices=['lsc', 'cpt', 'coj', 'elp', 'bpl'],
                        help="To which site to submit")
    parser.add_argument('--dome', default='doma', choices=['doma', 'domb', 'domc'], help="To which enclosure to submit")
    parser.add_argument('--telescope', default='1m0a')
    parser.add_argument('--instrument', default='fl12',
                        choices=['fl03', 'fl04', 'fl05', 'fl08', 'fl11', 'fl12', 'fl15', 'fl16', ],
                        help="To which instrumetn to submit")
    parser.add_argument('--start', default=None,
                        help="When to start x-talk calibration. If not given, defaults to \"NOW\"")
    parser.add_argument('--user', default='daniel_harbeck', help="Which user name to use for submission")
    parser.add_argument('--offsetRA', default=0, help="Extra pointing offset to apply R.A.")
    parser.add_argument('--offsetDec', default=0, help="Extra pointing offset to apply Dec")

    parser.add_argument('--CONFIRM', dest='opt_confirmed', action='store_true',
                        help='If set, block will be submitted.')

    parser.add_argument('--log_level', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
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
            _logger.error("Invalidt start time argument: ", args.start)
            exit(1)

    if ('auto' in args.name):
        # automatically find the best target
        args.name = getAutoCandidate(args)
        pass

    try:
        _logger.debug("Resolving target name")
        args.radec = SkyCoord.from_name(args.name)
    except:
        print("Resolving target name failed, giving up")
        exit(1)

    print("Resolved target %s at corodinates %s %s" % (args.name, args.radec.ra, args.radec.dec))
    return args


if __name__ == '__main__':
    args = parseCommandLine()

    createRequestsForStar(args)
