FROM ubuntu:bionic
LABEL maintainer="Nimbix, Inc."

# Update SERIAL_NUMBER to force rebuild of all layers (don't use cached layers)
ARG SERIAL_NUMBER
ENV SERIAL_NUMBER ${SERIAL_NUMBER:-20210113.1000}

ARG GIT_BRANCH
ENV GIT_BRANCH ${GIT_BRANCH:-master}

RUN apt-get -y update && \
    apt-get -y install curl && \
    curl -H 'Cache-Control: no-cache' \
        https://raw.githubusercontent.com/nimbix/image-common/$GIT_BRANCH/install-nimbix.sh \
        | bash -s -- --setup-nimbix-desktop --image-common-branch $GIT_BRANCH

RUN apt-get -y install default-jre

COPY nextflow/build/releases/nextflow-22.11.0-edge-all /opt/jarvice-nf/nextflow

WORKDIR /opt
COPY jarvice-nxf-submit-batch jarvice-nxf-qstatus jarvice.config jarvice-nf-run.sh /opt/jarvice-nf/
RUN ls -la /opt/jarvice-nf
ENV PATH="$PATH:/opt/jarvice-nf"

COPY NAE/screenshot.png /etc/NAE/screenshot.png
COPY NAE/AppDef.json /etc/NAE/AppDef.json
COPY NAE/help.html /etc/NAE/help.html

RUN curl --fail -X POST -d @/etc/NAE/AppDef.json https://cloud.nimbix.net/api/jarvice/validate

# Expose port 22 for local JARVICE emulation in docker
EXPOSE 22

# for standalone use
EXPOSE 5901
EXPOSE 443
