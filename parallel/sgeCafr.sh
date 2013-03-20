#!/bin/sh
#$ -j y
#$ -cwd
#$ -S /bin/bash
#$ -t 1-300
export INPUTFILE=$1
export JOBNAME=$2
Rscript parallel.sge.R $INPUTFILE $SGE_TASK_ID $SGE_TASK_LAST $JOBNAME
