=============================================
  cafr - parallel

		Wei-Yi Cheng
		2013-03-18
=============================================


The parallel directory contains all the scripts required for executing 
Attractor Finding algorithm in parallel. The files include:

sgeCafr.sh -- A shell script for running cafr in parallel under Sun Grid 
              Engine. Change the -t setting for the number of segments 
              the tasks are going to be divided.
parallel.sge.R -- An R script to run basic loading of datafiles, and 
                  perform parAttractorScanning based on the SGE setting.
mergeResults.R -- An R script to merge the results output by each worker.

An example to run cafr in parallel and then merge the result:

> qsub sgeCafr.sh <input_file_name> <job_name>

where <input_file_name> is the path to the input data file, which is an 
RData file contains a "ge" object with gene expression matrix. And 
<job_name> is an optional argument to name the output directory.

After the job has been finished, an output directory will be generated. 
To merge the output from each segment, run:

> Rscript mergeResults.R <output_directory>

Then it will generate two files:

attractorMatrix.rda -- an R object of a matrix with each row being an
                       attractor.
attractors.txt -- a tab-delimited txt file with every two columns being
                  the top genes of an attractor and their MI with the 
                  metagene.


