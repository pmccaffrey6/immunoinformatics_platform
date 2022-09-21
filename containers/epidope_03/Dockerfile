FROM nvcr.io/nvidia/tensorflow:19.11-tf1-py3
ENV VERSION 0.3
ENV TOOL epidope

RUN useradd -m epidopeuser

USER epidopeuser

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -P /home/epidopeuser && /bin/sh /home/epidopeuser/Miniconda3-latest-Linux-x86_64.sh -b -f

RUN /home/epidopeuser/miniconda3/bin/conda config --add channels conda-forge && \
    /home/epidopeuser/miniconda3/bin/conda config --add channels flomock && \
    /home/epidopeuser/miniconda3/bin/conda config --add channels default

RUN /home/epidopeuser/miniconda3/bin/conda install -y mamba
RUN /home/epidopeuser/miniconda3/bin/mamba install -y python=3.6 $TOOL=$VERSION
RUN /home/epidopeuser/miniconda3/bin/mamba clean -ay
RUN /home/epidopeuser/miniconda3/bin/conda install overrides==3.1 nltk==3.4

ENV PATH="${PATH}:/home/epidopeuser/miniconda3/bin"

ENTRYPOINT ["epidope"]
