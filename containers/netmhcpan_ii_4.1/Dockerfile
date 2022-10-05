FROM ubuntu:20.04

ENV PATH="${PATH}:/netMHCIIpan-4.1"

COPY container-build-items/netMHCIIpan-4.1a.Linux.tar.gz /netMHCIIpan-4.1a.tar.gz

RUN apt-get update && apt-get install -y wget tcsh gawk libenv-path-perl && \
    tar xzvf /netMHCIIpan-4.1a.tar.gz

COPY src/netMHCIIpan /netMHCIIpan-4.1/netMHCIIpan

RUN chmod +x /netMHCIIpan-4.1/netMHCIIpan
