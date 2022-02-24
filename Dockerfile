FROM python:3.10
MAINTAINER Las Cumbres Observatory <webmaster@lco.global>
# Configure application working directory
WORKDIR /lco/noisegainreport

RUN apt-get update -y \
        && apt-get install --no-install-recommends -y less vim\
        && apt-get clean -y \
        && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Python dependencies (libraries)
COPY requirements.txt .
RUN pip --no-cache-dir install --upgrade pip \
        && pip install --upgrade --force-reinstall setuptools\
        && pip --no-cache-dir install -r requirements.txt --upgrade
# Install application code
COPY . .
RUN python setup.py install
