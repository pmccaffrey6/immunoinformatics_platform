FROM ubuntu:20.04

COPY ../container-build-items/IEDB_BCell-3.1.tar.gz /
COPY ../container-build-items/netsurfp-2.0.Any.tar.gz /

RUN apt-get update && apt-get install nano curl tcsh gawk git python2 \
  && curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py \
  && python2 get-pip.py && ln -s /usr/bin/python2 /usr/bin/python \
  && apt-get install python3-dev python3-pip \
  && pip install -r bcell_standalone/requirements.txt
  && tar xzvf IEDB_BCell-3.1.tar.gz -C / && ./bcell_standalone/configure
  && mkdir /netsurfp && tar xzvf netsurfp-2.0.Any.tar.gz -C /netsurfp \
  && python3 /netsurfp/setup.py install

ENV NETSURFP_BIN=/usr/local/bin/netsurfp2

ENTRYPOINT ['python3', 'bcell_standalone/predict_antibody_epitope.py']
