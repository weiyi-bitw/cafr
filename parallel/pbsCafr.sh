#!/bin/bash

export INPUT=$1
export NTASK=$2
export JOBNAME=$3

echo "cd \$PBS_O_WORKDIR; Rscript parallel_pbs.R $INPUT $NTASK $JOBNAME" | 
  qsub -J 1-$NTASK -N $JOBNAME -q long -m a -j oe -M wei-yi.cheng@roche.com 
  #-l nodes=4
  #-l mem=$((NTASK * 2))gb,ncpus=$NTASK

