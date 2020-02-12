#!/usr/bin/env python3

'''
Simple Flask application to handle dynamically generating the HTML based on
the available images for the Long Term Photometric Zero Point web application
'''
import logging

from botocore.exceptions import ClientError
from flask import Flask
from flask import render_template

import collections
import datetime
import boto3
import os

app = Flask(__name__)

def get_last_update_time_from_object():
    client = boto3.client('s3')
    params = {
        'Bucket': os.environ.get('AWS_S3_BUCKET', None),
        'Key': 'gainhist-fa03full_frameNone.png',
    }
    try:
        response = client.head_object(**params)
        return response['LastModified']
    except:
        print ("Error while identifying last modified tag")
        pass
    return None



def s3_list_objects():
    '''
    A Python Generator which yields all Keys (filenames) in the
    configured S3 Bucket
    '''
    params = {
        'Bucket': os.environ.get('AWS_S3_BUCKET', None)
    }

    client = boto3.client('s3')
    paginator = client.get_paginator('list_objects_v2')
    for page in paginator.paginate(**params):
        for elem in page.get('Contents', []):
            yield elem



def build_site_info_dict():
    retval = collections.defaultdict(list)
    # TODO: instead of site -> filenames make a site->telescope->filenames dictionary
    for objdict in s3_list_objects():
        filename = objdict['Key']
        parts = filename.split('-')
        if len(parts) >= 2 and parts[0] == 'gainhist':
            cameracode = parts[1]
            retval[cameracode].append(filename)

    # Sort it appropriately
    for key, values in retval.items():
        values.sort(key = lambda x: x[0: 4])

    return retval

def generate_presigned_url(filename):
    # https://boto3.amazonaws.com/v1/documentation/api/latest/guide/s3-presigned-urls.html
    client = boto3.client('s3')
    params = {
        'Bucket': os.environ.get('AWS_S3_BUCKET', None),
        'Key': filename,
    }
    expiration = int(datetime.timedelta(days=7).total_seconds())
    try:
        response = client.generate_presigned_url('get_object', Params=params, ExpiresIn=expiration)
    except ClientError as ex:
        logging.error(str(ex))
        return None

    return response

@app.route('/')
@app.route('/index.html')
def indexhtml():
    params = {
        'generate_presigned_url': generate_presigned_url,
        'timestamp': get_last_update_time_from_object(),
        'info': build_site_info_dict(),
    }
    return render_template('index.html', **params)

@app.route('/healthz')
def healthz():
    return 'Healthy!\n'

def main():
    #app.run(port=8080)
    pass

if __name__ == '__main__':
    main()



# vim: set ts=4 sts=4 sw=4 et tw=120:
