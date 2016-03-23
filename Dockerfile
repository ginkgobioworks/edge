FROM ubuntu:14.04
MAINTAINER Giles Hall
RUN apt-get update && \
    apt-get install -y git build-essential libmysqlclient-dev \
        python python-dev python-setuptools python-virtualenv python-pip \
        nodejs nodejs-legacy npm primer3 ncbi-blast+
RUN npm install -g bower

# running environment
RUN useradd -ms /bin/bash edge
ADD / /home/edge/src/edge
RUN chown -R edge.edge /home/edge
USER edge
WORKDIR /home/edge/src/edge
RUN virtualenv /home/edge/env && \
    /home/edge/env/bin/python setup.py install && \
    /home/edge/env/bin/edge_manage.py assets build
EXPOSE 8000
ENTRYPOINT ["scripts/wait-for-it.sh", "-t", "60", "edgedb:3306", "--", "scripts/entrypoint.sh"]
CMD ["init", "run"]
