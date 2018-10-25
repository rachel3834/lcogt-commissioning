FROM python:3.6
MAINTAINER Las Cumbres Observatory <webmaster@lco.global>

EXPOSE 80

#RUN mkdir /home/archive  && /usr/sbin/groupadd -g 10000 "domainusers" \
#        && /usr/sbin/useradd -g 10000 -d /home/archive -M -N -u 10087 archive \
#        && chown -R archive:domainusers /home/archive



RUN apt-get update -y \
        && apt-get install --no-install-recommends -y cron supervisor less vim  parallel \
        && apt-get clean -y \
        && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt --upgrade

COPY . /lco/noisegain

COPY deploy/supervisor-app.conf /etc/supervisor/conf.d/
COPY deploy/crontab /etc/cron.d/
RUN chmod 0644 /etc/cron.d/crontab

CMD ["supervisord", "-n"]