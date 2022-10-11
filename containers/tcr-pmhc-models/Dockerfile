FROM continuumio/anaconda3:latest

ENV PATH="${PATH}:/usr/local/bin/scwrl4"
ENV KEY_MODELLER=MODELIRANJE

COPY container-build-items /container-build-items

RUN apt-get update && apt-get -y install hmmer wget gcc g++ && \
    cp -R /container-build-items/scwrl4 /usr/local/bin/scwrl4 && \
    cp /container-build-items/LYRA-1.0.ALL.tar.gz / && \
    tar xzvf /LYRA-1.0.ALL.tar.gz && \
    conda create -y -n TCRpMHCmodels python=2.7 && \
    /opt/conda/envs/TCRpMHCmodels/bin/pip install requests biopython==1.76 scikit-learn==0.15.2 scipy

RUN ( cd /LYRA && /opt/conda/envs/TCRpMHCmodels/bin/python /LYRA/setup.py install )

RUN conda install -n TCRpMHCmodels -c salilab modeller && \
    cp /container-build-items/TCRpMHCmodels-1.0.ALL.tar.gz / && \
    tar xzvf /TCRpMHCmodels-1.0.ALL.tar.gz

RUN ( cd /TCRpMHCmodels && /opt/conda/envs/TCRpMHCmodels/bin/python /TCRpMHCmodels/setup.py install )

#/opt/conda/envs/TCRpMHCmodels/bin/tcrpmhc_models /TCRpMHCmodels/tcrpmhc_models/test/TCRpMHC-test.fsa -n myname -p /tcrpmhc_output
