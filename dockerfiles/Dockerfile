# This is the barebones image for the GATK SV pipeline that includes only
# some basic bioinformatics utilities.

############### stage 0: build samtools and bcftools from source
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.9
ARG BCFTOOLS_VERSION=1.9
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
                 ca-certificates \
                 libbz2-dev \
                 libcurl4-openssl-dev \
                 liblzma-dev \
                 libncurses5-dev \
                 autoconf \
                 automake \
                 bzip2 \
                 gcc \
                 make \
                 wget \
                 zlib1g-dev && \
    wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && ./configure --without-curses --enable-libcurl && make -s all all-htslib && make install install-htslib && cd - && \
    rm -rf samtools-${SAMTOOLS_VERSION}* && \
    wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && ./configure --without-curses && make -s && make install && cd - && \
    rm -rf bcftools-${BCFTOOLS_VERSION}* && \
    apt-get -qqy purge autoconf automake bzip2 gcc make wget && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    samtools --help && \
    bcftools --help

############### stage 1: other commonly used bioinformatics utilities 
FROM ubuntu:20.04

ENV HOME=/root

# Define default command.
CMD ["bash"]

# copy from previous stage the binaries from samtools build
COPY --from=0 /usr/local/bin/* /usr/local/bin/

ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install --no-install-recommends \
        apt-utils \
        bedtools \
        build-essential \
        ca-certificates \
        curl \
        gawk \
        git \
        less \
        libpcre3-dev \
        liblzma-dev \
        libbz2-dev \
        libicu-dev \
        openjdk-8-jdk \
        python3-dev \
        python3-pip \
        python3-setuptools \
        r-base \
        r-base-dev \
        tabix \
        vcftools \
        wget \
        xvfb \
        zip \
        zlib1g-dev

# Install Google Cloud SDK
RUN apt-get update && apt-get install -y apt-transport-https ca-certificates gnupg && \
    pip3 install -q -U --no-cache-dir crcmod wheel && \
    echo "deb [signed-by=/usr/share/keyrings/google-cloud.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -s https://packages.cloud.google.com/apt/doc/apt-key.gpg | gpg --dearmor > /usr/share/keyrings/google-cloud.gpg && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends google-cloud-sdk && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image

# Clean up
RUN apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_*

# install vapor
RUN pip3 install --no-cache-dir Cython && \
    cd /opt && \
    git clone https://github.com/kjaisingh/vapor.git && \
    cd /opt/vapor/ && \
    git checkout fe8ee37adcfc2fecd84b6f8a83933edacc1901bf && \
    pip3 install -e .
