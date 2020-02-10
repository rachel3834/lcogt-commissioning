FROM python:3.7

# Configure application working directory
WORKDIR /lco/noisegainreport

# Install Python dependencies (libraries)
COPY requirements.txt .
RUN pip --no-cache-dir install -r requirements.txt --upgrade

# Install application code
COPY . .
RUN python setup.py install
