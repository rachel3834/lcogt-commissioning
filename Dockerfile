FROM python:3.7

# Configure application working directory
WORKDIR /lco/noisegainreport

RUN apt-get update -y \
        && apt-get install --no-install-recommends -y less vim postgresql-common postgresql-server-dev-11\
        && apt-get clean -y \
        && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install Python dependencies (libraries)
COPY requirements.txt .
RUN pip --no-cache-dir install -r requirements.txt --upgrade

# Install application code
COPY . .
RUN python setup.py install
