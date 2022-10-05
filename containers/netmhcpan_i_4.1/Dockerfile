FROM ubuntu:20.04

ENV PATH="${PATH}:/netMHCpan-4.1"

COPY container-build-items/netMHCpan-4.1b.Linux.tar.gz /netMHCpan-4.1b.tar.gz

RUN apt-get update && apt-get install -y wget tcsh gawk && \
    tar xzvf /netMHCpan-4.1b.tar.gz && \
    wget https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz && \
    tar xzvf /data.tar.gz

COPY src/netMHCpan /netMHCpan-4.1/netMHCpan

RUN chmod +x /netMHCpan-4.1/netMHCpan && \
    rm /netMHCpan-4.1/Linux_x86_64/data && \
    ln -s /data /netMHCpan-4.1/Linux_x86_64/data
