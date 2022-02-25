import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lcocommissioning",
    version="2.1.0",
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
                            'submit_muscat_observation = lcocommissioning.muscat.submitMuscatObservation:main',
                            'submit_cdk_observation = lcocommissioning.cdk.submit_cdk_observation:main',
                            'get_sbigframe = lcocommissioning.sbig.get_sbig_frame:main'],

    },
    zip_safe=False,

)
