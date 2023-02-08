#!/bin/bash
export JARVICE_API_URL=$2
export JARVICE_USER=$4
export JARVICE_API_KEY=$6
export NXF_HOME=/data/.nextflow

cd /data

echo -e "\e[1;34m----------------------------------\nWelcome to JARVICE Nextflow Shell\n----------------------------------\e[0m\nUse \e[1;31mjarvice-nf-run.sh <pipeline>\e[0m to launch\n"
/bin/bash -il
