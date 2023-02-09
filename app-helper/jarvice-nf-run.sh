#!/bin/bash
export NXF_HOME=/tmp/.nextflow
nextflow run $1 -c /opt/jarvice-nf/jarvice.config
