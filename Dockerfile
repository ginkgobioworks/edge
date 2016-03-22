FROM ubuntu:14.04
MAINTAINER Giles Hall
RUN apt-get update && \
    apt-get install -y python python-dev build-essential python-setuptools python-virtualenv

# running environment
ENV EDGE_DEFAULT_DB=sqlite
ENV EDGE_TESTING=1
ADD / /opt/edge
WORKDIR /opt/edge
RUN virtualenv env && ./env/bin/python setup.py install
ENTRYPOINT ["/opt/edge/startup.sh"]

# MySQL
#VOLUME ["/var/lib/mysql"]
#RUN apt-get install -y mysql-server
