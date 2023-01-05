FROM ubuntu:xenial
LABEL maintainer="Nimbix, Inc."

# Update SERIAL_NUMBER to force rebuild of all layers (don't use cached layers)
ARG SERIAL_NUMBER
ENV SERIAL_NUMBER ${SERIAL_NUMBER:-20180124.1405}

ARG GIT_BRANCH
ENV GIT_BRANCH ${GIT_BRANCH:-main}

RUN apt-get -y update && apt-get -y install curl

RUN apt-get -y install ncbi-blast+


# Make a dummy /data dir so the symlink doesnt dangle
RUN mkdir /data

# Point the vault dir to the /nextflow dir
RUN ln -s /data /nextflow && chmod 777 /nextflow


# Expose port 22 for local JARVICE emulation in docker
EXPOSE 22

# for standalone use
EXPOSE 5901
EXPOSE 443
