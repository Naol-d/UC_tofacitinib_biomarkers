FROM ubuntu:22.04

ENV CELLRANGER_VERSION=9.0.1
ENV PATH=/opt/cellranger-${CELLRANGER_VERSION}:$PATH

RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    build-essential \
    dpkg-dev \
    gnupg \
    procps \
    python3 \
    python3-pip \
    python3-venv \
    locales \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*


# NOTE: wget link expires on the 17th August 2025, removed for security purposes
RUN wget -O /opt/cellranger-${CELLRANGER_VERSION}.tar.gz \
    "<ENTER CELLRANGER wget link>" \
 && cd /opt && \
    tar -xzf cellranger-${CELLRANGER_VERSION}.tar.gz && \
    rm cellranger-${CELLRANGER_VERSION}.tar.gz

