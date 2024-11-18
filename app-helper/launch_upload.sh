#!/bin/bash

while [[ $# -gt 0 ]]; do
    case "$1" in
        -APIURL)
            JARVICE_API_URL="$2"
            shift 2
            ;;
        -APIUSER)
            JARVICE_USER="$2"
            shift 2
            ;;
        -APIKEY)
            JARVICE_API_KEY="$2"
            shift 2
            ;;
        -VAULT)
            JARVICE_VAULT="$2"
            shift 2
            ;;
        *)
            shift 2
            ;;
    esac
done

export JARVICE_API_URL
export JARVICE_USER
export JARVICE_API_KEY
export JARVICE_VAULT

echo "JARVICE_API_URL=$JARVICE_API_URL"
echo "JARVICE_USER=$JARVICE_USER"
echo "JARVICE_API_KEY=$JARVICE_API_KEY"
echo "JARVICE_VAULT=$JARVICE_VAULT"

export NXF_HOME=/tmp/.nextflow

cd /data

echo -e "\e[1;34m----------------------------------\nWelcome to JARVICE Nextflow Shell\n----------------------------------\e[0m\n"

echo " Starting Nextflow execution..."

jarvice-nf-run.sh /opt/flow.nf
