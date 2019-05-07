import os
import requests

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


# Tool to submit a user request to LCO Scheduler via Valhalla interface
def send_to_scheduler(user_request, dosubmit=False):
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