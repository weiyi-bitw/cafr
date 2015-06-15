LoadTextMat <- function(file, sep='\t'){
  line <- readLines(file)
  tokens <- strsplit(line[2], sep)[[1]]
  n <- length(tokens)-1
  m <- length(line)-1
  mset <- matrix(NA, m, n)
  tokens <- strsplit(line[1], sep)[[1]]
  if(length(tokens)==n){
    cname <- tokens[1:n]
  }else if(length(tokens)==n+1){
    cname <- tokens[2:(n+1)]
  }else{
    stop("The number of entries in header is incompatible with the number of other lines!")
  }
  rname <- rep(NA, m)
  b <- txtProgressBar(style=3)
  for(i in 1:m){
    tokens <- strsplit(line[i+1], sep)[[1]]
    tokens[tokens=="[Not Available]" | tokens=="null" | tokens=="NA"] <- NA
    rname[i] <- tokens[1]
    mset[i,] <- tokens[2:(n+1)]
    if(i %% 100 == 0){
      setTxtProgressBar(b, i/m)
    }
  }
  cat("\n")
  colnames(mset) <- cname
  rownames(mset) <- rname
  mset
}
