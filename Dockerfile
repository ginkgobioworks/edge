# Development Dockerfile to make testing easier under a standardized environment.
# XXX Do not use for production as-is

FROM python:3.6-buster
LABEL maintainer Ginkgo Bioworks <devs@ginkgobioworks.com>

ARG GIT_USER_NAME="Curious Default"
ARG GIT_USER_EMAIL="devs@ginkgobioworks.com"

RUN git config --global user.name "$GIT_USER_NAME" \
    && git config --global user.email "$GIT_USER_EMAIL"


ARG DEBIAN_FRONTEND=noninteractive
RUN curl -sL https://deb.nodesource.com/setup_12.x | bash -
RUN apt-get update && apt-get install --assume-yes --verbose-versions \
  apt-utils \
  default-mysql-client \
  nodejs \
  ncbi-blast+ \
  primer3

RUN npm install --global bower
RUN echo '{ "allow_root": true }' > $HOME/.bowerrc

ENV EDGE_HOME=/usr/src/edge
RUN mkdir -p $EDGE_HOME
WORKDIR $EDGE_HOME

COPY requirements-core.txt ./
RUN pip install -r requirements-core.txt

COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY . ./
RUN pip install --editable .

RUN mkdir -p $EDGE_HOME/ncbi/blastdb
ENV EDGE_IP=0.0.0.0 \
    EDGE_PORT=8000 \
    NCBI_BIN_DIR=/usr/bin \
    PRIMER3_BIN=/usr/bin/primer3_core \
    PRIMER3_CONFIG_DIR=/etc/primer3_config/

EXPOSE $EDGE_PORT
CMD ["make", "start"]
