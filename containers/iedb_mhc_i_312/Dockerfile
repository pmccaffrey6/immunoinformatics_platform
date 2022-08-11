FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

COPY container-build-items/IEDB_MHC_I-3.1.2.tar.gz /

RUN apt-get update && apt-get install -y nano gawk tcsh python3-dev python3-pip \
  && tar xzvf IEDB_MHC_I-3.1.2.tar.gz -C / \
  && ln -s /usr/bin/python3 /usr/bin/python

WORKDIR /mhc_i

RUN ./configure
