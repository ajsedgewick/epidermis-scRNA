FROM ubuntu:16.04
RUN apt-get update \
    && apt-get install -y python3 python3-pip libreadline-dev python3-tk less git \
    && pip3 install numpy scipy pandas sklearn argparse \
    && git clone git://github.com/pkathail/magic.git \
    && cd magic \
    && git checkout 3f5db2c193699c32974bb6f738149ad0b1fcef0b \
    && pip3 install . \
    && rm -rf /var/lib/apt/lists/*	

ENV R_BASE_VERSION 3.4.2*
ENV DISTRIBUTION_CODENAME xenial

RUN  apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E084DAB9 \ 
    && echo "deb http://cran.r-project.org/bin/linux/ubuntu ${DISTRIBUTION_CODENAME}/" >> /etc/apt/sources.list

RUN apt-get update \
    && apt-get install -y --no-install-recommends r-base-core=${R_BASE_VERSION}* \
               r-base=${R_BASE_VERSION}* \
       	       r-base-dev=${R_BASE_VERSION}* \
               r-recommended=${R_BASE_VERSION}* \
               r-cran-kernsmooth r-cran-mgcv=1.8-22* r-cran-mass=7.3-47* \
               r-cran-class \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update \ 
    && apt-get install -y --no-install-recommends \
       python-software-properties \
       ed \
       less \
       locales \
       vim-tiny \
       wget \
       ca-certificates \
       software-properties-common \
       libssl-dev \
       libcurl4-openssl-dev \
       libxml2-dev \
       libgsl0-dev \
       libcairo2-dev \
       libhdf5-dev \
   && apt-get update 

# RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
#    && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default


RUN Rscript -e "install.packages(c('devtools'), repos='http://cran.us.r-project.org')"

RUN Rscript -e "install.packages(c('Matrix', 'ggplot2', 'dplyr', 'ggvis', 'feather', 'openxlsx', 'Seurat', 'statmod', 'svglite', 'reshape'), repos='http://cran.us.r-project.org')"

RUN Rscript -e "source('https://bioconductor.org/biocLite.R'); \
                biocLite('limma'); \
                biocLite('zinbwave'); \
                biocLite('clusterExperiment')"

RUN Rscript -e "library(devtools); install_github('igraph/rigraph'); \
                install_github('kstreet13/slingshot')"


RUN pip3 install feather-format

# can specify UID here to make your life easier
RUN useradd scuser && usermod -a -G sudo scuser

RUN mkdir -p  /home/scuser \
     && echo "R_MAX_NUM_DLLS=500" > /home/scuser/.Renviron \
     && chown -R scuser /home/scuser \
     && chgrp -R scuser /home/scuser

USER scuser
