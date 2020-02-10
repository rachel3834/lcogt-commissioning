FROM python:3.7
MAINTAINER Las Cumbres Observatory <webmaster@lco.global>


#RUN mkdir /home/archive  && /usr/sbin/groupadd -g 10000 "domainusers" \
#        && /usr/sbin/useradd -g 10000 -d /home/archive -M -N -u 10087 archive \
#        && chown -R archive:domainusers /home/archive

WORKDIR /lco/noisegainreport

RUN apt-get update -y \
        && apt-get install --no-install-recommends -y cron less vim \
        && apt-get clean -y \
        && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


COPY requirements.txt .
RUN pip --no-cache-dir install --upgrade pip \
        && pip --no-cache-dir install -r requirements.txt --upgrade

COPY . .
RUN python setup.py install


#COPY deploy/supervisor-app.conf /etc/supervisor/conf.d/
#COPY deploy/crontab /etc/cron.d/
#RUN chmod 0644 /etc/cron.d/crontab

#CMD ["supervisord", "-n"]