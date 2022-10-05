FROM ubuntu:20.04

ENV GALAXY_HOME=/GalaxyPepDock-master
ENV EXEC_MPI=/hydra-3.2/mpiexec.hydra
ENV NSLOTS=10

COPY container-build-items/GalaxyPepDock.zip /GalaxyPepDock.zip

RUN apt-get update && apt-get install -y unzip zip wget nano gcc g++ make && \
    unzip /GalaxyPepDock.zip && \
    chmod +x /GalaxyPepDock-master/install.sh && \
    ./GalaxyPepDock-master/install.sh && \
    chmod +x ./GalaxyPepDock-master/bin/GalaxyPepDock.ubuntu1604 && \
    wget http://www.mpich.org/static/downloads/3.2/hydra-3.2.tar.gz && \
    tar xzvf hydra-3.2.tar.gz && \
    ./hydra-3.2/configure &&
