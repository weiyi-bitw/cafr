#!/bin/bash

export INPUT=$1
export JOBID=$2

echo "cd \$PBS_O_WORKDIR; Rscript pbs_geupgma.R $INPUT $JOBID" | 
  qsub -N $JOBID -q long -m a -j oe -M wei-yi.cheng@roche.com -l mem=12gb

