# Development Dockerfile to make testing easier under a standardized environment.
# XXX Do not use for production as-is

FROM python:2.7
LABEL maintainer Ginkgo Bioworks <devs@ginkgobioworks.com>

ARG GIT_USER_NAME="Curious Default"
ARG GIT_USER_EMAIL="devs@ginkgobioworks.com"

RUN git config --global user.name "$GIT_USER_NAME" \
    && git config --global user.email "$GIT_USER_EMAIL"

ENV INSTALL_ROOT=/usr/src

WORKDIR $INSTALL_ROOT

ARG install_ncbi=true
COPY ncbi ncbi
ENV NCBI_HOME=$INSTALL_ROOT/ncbi
RUN if $install_ncbi; then cd $NCBI_HOME && ./install; else echo "Skipping NCBI BLAST install"; fi

ARG install_primer3=true
COPY primer3 primer3
ENV PRIMER3_HOME=$INSTALL_ROOT/primer3
RUN if $install_primer3; then cd $PRIMER3_HOME && ./install; else echo "Skipping Primer3 install"; fi


RUN apt-get update
RUN apt-get install --assume-yes apt-utils nodejs nodejs-legacy npm

RUN npm install --global bower
RUN echo '{ "allow_root": true }' > $HOME/.bowerrc


ENV EDGE_HOME=$INSTALL_ROOT/edge
RUN mkdir -p $EDGE_HOME
WORKDIR $EDGE_HOME

COPY requirements-core.txt ./
RUN pip install -r requirements-core.txt

COPY requirements.txt ./
RUN pip install -r requirements.txt

COPY . ./
RUN pip install --editable .

RUN \
  if $install_ncbi; then \
    mkdir -p $EDGE_HOME/ncbi/blastdb; \
    cd $EDGE_HOME/ncbi; \
    ln -sf $NCBI_HOME/ncbi-blast-2.2.27+-src/c++/GCC492-Debug64/bin bin; \
  fi

RUN \
  if $install_primer3; then \
    mkdir -p $EDGE_HOME/primer3; \
    cd $EDGE_HOME/primer3; \
    ln -sf $PRIMER3_HOME/primer3_core .; \
    ln -sf $PRIMER3_HOME/primer3_config .; \
  fi

ENV EDGE_IP=0.0.0.0 \
    EDGE_PORT=8000
EXPOSE $EDGE_PORT
CMD ["make", "start"]
