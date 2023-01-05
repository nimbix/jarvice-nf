#!/bin/bash

sshfs -o IdentityFile=/home/rep/.ssh/google_compute_engine rep@sandbox.us-central1-a.product-engineering-devtest:/data/rep /nextflow
cp -v /home/rep/projects/nimbix/jarvice-nf/* /nextflow/

ls -la /nextflow
