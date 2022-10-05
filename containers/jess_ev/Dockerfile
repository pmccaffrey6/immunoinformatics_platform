FROM continuumio/anaconda3

COPY src /

ENV PATH="${PATH}:/gurobi952/linux64/bin"
ENV WLSACCESSID=931c9283-46b4-4af0-b965-2620305e3695
ENV WLSSECRET=d6e250d9-da06-449c-b69e-fa66089a76ef
ENV LICENSEID=876137
ENV GRB_LICENSE_FILE=/gurobi.lic
ENV GUROBI_HOME=/gurobi952/linux64/

RUN wget https://packages.gurobi.com/9.5/gurobi9.5.2_linux64.tar.gz && \
    tar xzvf gurobi9.5.2_linux64.tar.gz && \
    conda create -y -n jessev --file /JessEV/packages.txt -c conda-forge
