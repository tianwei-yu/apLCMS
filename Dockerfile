FROM ubuntu:20.04

WORKDIR /usr/src/recetox-aplcms
COPY conda/environment-dev.yaml /usr/src/recetox-aplcms/conda/environment-dev.yaml

# install R-base dependencies and wget
RUN apt-get update -y && \
    apt-get install -y libxml2-dev && \
    apt-get install -y libssl-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libcgal-dev && \
    apt-get install -y libglu1-mesa-dev && \
    apt-get install -y wget

# download and install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /bin/local/miniconda
ENV PATH /bin/local/miniconda/bin:$PATH

# configure conda and create environment
RUN conda init
RUN conda update conda
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict && \
    conda env create -f conda/environment-dev.yaml

ENTRYPOINT ["/bin/local/miniconda/envs/recetox-aplcms/bin/R", "-e", "devtools::test()"]