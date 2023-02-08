#!/bin/bash
export NXF_HOME=/data/.nextflow
nextflow run $1 -c /opt/jarvice-nf/jarvice.config
