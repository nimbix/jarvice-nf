#!/bin/bash
umask 000
sudo chown -R nimbix:nimbix /data
nextflow run $1 -c /opt/jarvice-nf/jarvice.config
