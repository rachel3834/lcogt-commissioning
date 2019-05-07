import logging
import os
import math
import ephem
import requests
import datetime as dt

_log = logging.getLogger(__name__)
# LCO Request submisison definitions
LAKE_URL = 'http://lake.lco.gtn'
VALHALLA_URL = os.getenv('VALHALLA_URL', 'http://valhalla.lco.gtn')
VALHALLA_API_TOKEN = os.getenv('VALHALLA_API_TOKEN', '')

# LCO sites
lco_site_lonlat = {'bpl': (-119.863103, 34.433161),
                   'coj': (149.0708466, -31.2728196),
                   'cpt': (20.8124, -32.3826),
                   'elp': (-104.015173, 30.679833),
                   'lsc': (-70.8049, -30.1673666667),
                   'ogg': (-156.2589, 34.433161),
                   'sqa': (-120.04222167, 34.691453333),
                   'tfn': (-16.511544, 28.300433), }

# Dictionary of NRES instances
nres_instruments = {'lsc': 'nres01',
                    'elp': 'nres02',
                    'cpt': 'nres03',
                    'tlv': 'nres04',
                    }

def getEphemObForSiteAndTime (sitecode, dateobs):
    site = ephem.Observer()
    lon, lat = lco_site_lonlat[sitecode]
    site.lat = lat * math.pi / 180
    site.lon = lon * math.pi / 180
    site.date = ephem.Date(dateobs)
    return site

def is_valid_lco_site (sitecode):
    return (sitecode in lco_site_lonlat)

def send_to_scheduler(user_request, dosubmit=False):
    """Submit a user request to LCO Scheduler via Valhalla interface
    """

    auth = 'Token {token}'.format(token=VALHALLA_API_TOKEN)
    print(auth)
    url = '{api_root}/api/userrequests/'.format(api_root=VALHALLA_URL)
    headers = {
        'Authorization': auth
    }
    validation_check = requests.post(url + 'validate/', json=user_request, headers=headers).json()
    print(validation_check)
    if not validation_check['errors']:
        print('There are no errors.')
        if dosubmit:
            submitted = requests.post(url, json=user_request, headers=headers).json()
            print('Submitted request information: {}'.format(submitted))
    else:
        print('Output of the validation check: {}'.format(validation_check))
        print('UserRequest that was validated: {}'.format(user_request))


def send_to_lake(block, dosubmit=False):
    """ Submit a user request to LCO POND"""
    if dosubmit:
        response = requests.post(LAKE_URL + '/blocks/', json=block)
        try:
            response.raise_for_status()
            _log.info(
                'Submitted block with id: {0}. Check it at {1}/blocks/{0}'.format(response.json()['id'], LAKE_URL))
        except Exception:
            _log.error(
                'Failed to submit block: error code {}: {}'.format(response.status_code, response.content))

