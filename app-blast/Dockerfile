# This is a generic container with NCBI blast installed, to use with Jarvice /batch endpoint
# TODO - make this slim or reuse existing public image

FROM ubuntu:xenial
LABEL maintainer="Nimbix, Inc."

# Update SERIAL_NUMBER to force rebuild of all layers (don't use cached layers)
ARG SERIAL_NUMBER
ENV SERIAL_NUMBER ${SERIAL_NUMBER:-20180124.1405}

ARG GIT_BRANCH
ENV GIT_BRANCH ${GIT_BRANCH:-main}

RUN apt-get -y update && apt-get -y install curl

RUN apt-get -y install ncbi-blast+

# Expose port 22 for local JARVICE emulation in docker
EXPOSE 22

# for standalone use
EXPOSE 5901
EXPOSE 443
