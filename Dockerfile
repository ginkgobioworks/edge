FROM python:2.7

WORKDIR /
COPY ncbi ncbi
WORKDIR /ncbi
RUN ./install

WORKDIR /
COPY primer3 primer3
WORKDIR /primer3
RUN ./install

WORKDIR /
COPY requirements.txt requirements.txt
COPY requirements-dev.txt requirements-dev.txt
RUN pip install -r requirements-dev.txt

WORKDIR /
ENV EDGE_HOME /edge
RUN mkdir -p $EDGE_HOME
COPY . $EDGE_HOME

RUN mkdir -p $EDGE_HOME/ncbi/blastdb
WORKDIR $EDGE_HOME/ncbi
RUN ln -sf /ncbi/ncbi-blast-2.2.27+-src/c++/GCC492-Debug64/bin bin

WORKDIR $EDGE_HOME/primer3
RUN ln -sf /primer3/primer3_core .
RUN ln -sf /primer3/primer3_config .

WORKDIR $EDGE_HOME/src
CMD python manage.py runserver 0.0.0.0:8000
