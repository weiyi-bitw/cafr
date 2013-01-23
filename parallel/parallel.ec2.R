# R script for Amazon EC2 cluster
#
# arguments for script:
#	args[1] : input gene expression data matrix
#
#

require(parallel)
require(caret)

# Install cafr package
installCAFR = function(){
	download.file("http://dl.dropbox.com/u/11986954/cafr/cafr_0.22.tar.gz", destfile="./cafr_0.22.tar.gz")
	install.packages("./cafr_0.22.tar.gz", repos=NULL)
	library(cafr)
	unlink("./cafr_0.22.tar.gz")
}

args <- commandArgs(TRUE)

#=================================
# Create socket cluster 
#
# code modified from Bioconductor cloud AMI tutorial:
#	http://www.bioconductor.org/help/bioconductor-cloud-ami/
#
#=================================
cat("Cluster initialization...\n");flush.console()

lines <- readLines("/usr/local/Rmpi/hostfile.plain")
hosts <- character()
for (line in lines)
{
    x <- (strsplit(line[[1]], " "))
    hosts <-
        c(hosts, rep.int(x[[1]][1], as.integer(x[[1]][2])))
}

cl <- makePSOCKcluster(hosts,
	master=system("hostname -i", intern=TRUE))

numWorkers = length(hosts)-1
widList = list()
for(i in 1:numWorkers){
	widList = c(widList, i)
}
#================================
# Load data
#===============================

ge = load.exp(args[1])
m = nrow(ge)
taskList = createFolds(1:m, k=numWorkers)

attractorOut = parLapply(cl, widList, 
  function(ge, tl, wid, a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, verbose=FALSE){
	parAttractorScanning(ge, tl, wid, a, maxIter, epsilon, bin, so)
  },
  ge=ge,
  tl=taskList
)


