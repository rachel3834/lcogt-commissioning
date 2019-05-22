import lcocommissioning.common.common as common
import datetime as dt
import logging
_logger = logging.getLogger(__name__)
logging.basicConfig(level=getattr(logging, 'DEBUG'),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
start = dt.datetime.utcnow() + dt.timedelta(minutes=1)
end = start + dt.timedelta(minutes=2)
start = str(start).replace(' ', 'T')
end = str(end).replace(' ', 'T')


request = {

    "molecules": [
        {

            "tag_id": "LCOGT",
            "user_id": "daniel_harbeck",
            "prop_id": "engineering test",
            "group": "testbias",
            "exposure_count": 1,
            "exposure_time": 30,
            "bin_x": 1,
            "bin_y": 1,
            "inst_name": "nres04",
            "priority": 1,
            "type": "NRES_DARK",
            "ag_filter": "",
            "exposure_time": "60",

        }
    ],
    "start": start,
    "end": end,
    "site": "tlv",
    "observatory": "igla",
    "telescope": "1m0a",
    "instrument_class": "1M0-NRES-SCICAM",
    "length": 0,
    "priority": 10,


}

common.send_to_lake (request, dosubmit=True)