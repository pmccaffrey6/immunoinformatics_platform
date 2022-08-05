FROM ubuntu:20.04

COPY /discotope-1.1 /discotope-1.1

RUN apt-get update && apt-get install -y nano python python-dev python3 curl gcc \
  && mkdir /usr/local/python && mkdir /usr/local/python/bin \
  && ln -s /usr/bin/python /usr/local/python/bin/python \
  && curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py \
  && python get-pip.py && pip install numpy && pip install biopython==1.55

ENTRYPOINT ["/discotope-1.1/discotope"]
