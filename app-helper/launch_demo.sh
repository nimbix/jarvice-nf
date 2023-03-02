#!/bin/bash
export JARVICE_API_URL=$2
export JARVICE_USER=$4
export JARVICE_API_KEY=$6
export JARVICE_VAULT=$8
export JARVICE_MACHINETYPE=${10}
export NXF_HOME=/tmp/.nextflow

export PATH=/opt/jarvice-nf/:$PATH

echo -e "\e[1;34m----------------------------------\nWelcome to JARVICE Nextflow Shell\n----------------------------------\e[0m\n"

echo " Starting Nextflow execution..."

mkdir -p /data/blast
cd /data/blast
$(which cp) -f /opt/samples/* .
sed -i "s/n1/$JARVICE_MACHINETYPE/" jarvice_blast.nf
echo "Submitting demo blast job with following content:"
echo "-------------------------------------------------"
cat jarvice_blast.nf
echo "-------------------------------------------------"
/opt/jarvice-nf/jarvice-nf-run.sh jarvice_blast.nf
