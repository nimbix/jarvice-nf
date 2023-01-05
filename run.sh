#!/bin/bash

export PATH=$PATH:/data/bin
./nextflow/launch.sh run $1 -c /nextflow/jarvice.config
