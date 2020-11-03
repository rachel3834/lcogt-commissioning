import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lcocommissioning",
    version="2.0.14",
    author="Daniel Harbeck",
    author_email="dharbeck@lco.global",
    description="Tool to characterize CCD detectors and other commissioning tasks for the LCO observatory.",
    long_description=long_description,
    url="https://github.com/LCOGT/lcogt-commissioning",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    entry_points = {
        'console_scripts': ['submitNRESObservation = lcocommissioning.nres.submitNRESObservation:main',
                            'noisegainmef = lcocommissioning.noisegainrawmef:main',
                            'submitXtalkObservation = lcocommissioning.submitXtalkObservation:main',
                            'analysegainhistory = noisegaincrawler.analysegainhistory:main',
                            'crawlnoisegain = noisegaincrawler.crawl_noisegain:main',
                            'submitNamedModeTest = lcocommissioning.submitNameModeTest:main',
                            'sinistrocrosstalk = lcocommissioning.crosstalk:main',
                            'focuscurve = lcocommissioning.focuscurve:main',
                            'effocus = lcocommissioning.ef_focuscalibration:main',
                            'submit_floyds_calibration = lcocommissioning.floyds.submitFloydsCalibration:main',
                            'submit_floyds_observation = lcocommissioning.floyds.submitFloydsObservation:main',
                            'submit_muscat_observation = lcocommissioning.muscat.submitMuscatObservation:main'
                            ],

    },
    install_requires=[
        "astropy==4.0",
        "numpy==1.18.1",
        "matplotlib==3.3.1",
        "ephem==3.7.7.0",
        "requests==2.22.0",
        "scipy==1.4.1",
        "sep==1.0.3",
        "peakutils==1.3.3",
        "sqlalchemy==1.3.13",
        "sqlalchemy_utils==0.36.1",
        "elasticsearch==7.5.1",
        "elasticsearch_dsl==7.1.0",
        "boto3==1.11.14",
        "Flask==1.1.1",
        "gunicorn[gevent]==20.0.4",
        "psycopg2-binary==2.8.4",
    ],
)
