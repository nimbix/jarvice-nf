#!/bin/bash
export JARVICE_API_URL=$2
export JARVICE_USER=$4
export JARVICE_API_KEY=$6
export JARVICE_VAULT=$8
export NXF_HOME=/tmp/.nextflow

cd /data

export PATH=/opt/jarvice-nf:$PATH

echo -e "\e[1;34m----------------------------------\nWelcome to JARVICE Nextflow Shell\n----------------------------------\e[0m\n"

echo " Starting Nextflow execution..."

/opt/jarvice-nf/jarvice-nf-run.sh /opt/flow.nf
