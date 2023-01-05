#!/bin/bash
umask 000
nextflow run $1 -c /opt/jarvice-nf-master/jarvice.config
