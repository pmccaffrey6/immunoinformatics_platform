FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

COPY container-build-items/IEDB_BCell-3.1.tar.gz /
COPY container-build-items/netsurfp-2.0.Any.tar.gz /

RUN apt-get update && apt-get install -y wget nano curl tcsh gawk git python2 software-properties-common pkg-config apt-utils libfreetype6-dev \
  && curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py \
  && curl https://bootstrap.pypa.io/pip/3.6/get-pip.py --output get-pip3.py \
  && add-apt-repository -y ppa:deadsnakes/ppa && apt-get update && apt-get -y install python3.6 python3.6-distutils \
  && python2 get-pip.py && ln -s /usr/bin/python3.6 /usr/bin/python \
  && apt-get install -y python3.6-dev && python3.6 /get-pip3.py \
  && python3.6 -m pip install pip==19.3.1 \
  && tar xzvf /IEDB_BCell-3.1.tar.gz -C /

WORKDIR /bcell_standalone

RUN pip install -r requirements.txt && /bin/sh configure

RUN mkdir /netsurfp && tar xzvf /netsurfp-2.0.Any.tar.gz -C /netsurfp
WORKDIR /netsurfp
RUN pip3.6 install numpy==1.19.2 scipy sklearn && python3.6 /netsurfp/setup.py install

ENV NETSURFP_BIN=/usr/local/bin/netsurfp2

ENTRYPOINT ["python", "/bcell_standalone/predict_antibody_epitope.py"]
