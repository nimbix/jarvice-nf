FROM ubuntu:focal
LABEL maintainer="Nimbix, Inc."

ARG GIT_BRANCH
ENV GIT_BRANCH ${GIT_BRANCH:-main}

ENV JARVICE_DESKTOP_ITER 1
RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get -y install ca-certificates curl --no-install-recommends

RUN curl -H 'Cache-Control: no-cache' \
    https://raw.githubusercontent.com/nimbix/jarvice-desktop/master/install-nimbix.sh \
    | bash -s

RUN DEBIAN_FRONTEND=noninteractive apt-get -y install default-jre

RUN java -version

COPY nextflow/build/releases/nextflow-22.11.0-edge-all /opt/jarvice-nf/nextflow

WORKDIR /opt
COPY jarvice-nxf-submit-batch jarvice-nxf-qstatus jarvice.config jarvice-nf-run.sh launch.sh /opt/jarvice-nf/
RUN ls -la /opt/jarvice-nf
ENV PATH="$PATH:/opt/jarvice-nf"

COPY NAE/screenshot.png /etc/NAE/screenshot.png
COPY NAE/AppDef.json /etc/NAE/AppDef.json

RUN curl --fail -X POST -d @/etc/NAE/AppDef.json https://cloud.nimbix.net/api/jarvice/validate
