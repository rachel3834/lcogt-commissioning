import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lcocommissioning",
    version="1.0.1",
    author="Daniel Harbeck",
    author_email="dharbeck@lco.global",
    description="Tool to characterize CCD detectors and other commissioning tasks for the LCO observatory.",
    long_description=long_description,
    url="https://github.com/drhaz/lcogt-commissioning",
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
                            'analysgainhistory = lcocommissioning.analysegainhistory:main',
                            'crawlnoisegain = lcocommissioning.crawl_noisegain:main',
                            'submitNamedModeTest = lcocommissioning.submitNameModeTest:main',
                            'sinistrocrosstalk = lcocommissioning.crosstalk:main',
                            'focuscurve = lcocommissioning.focuscurve:main'],

    }
)