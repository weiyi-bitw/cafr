ExtractGeneNameTCGA <- function(identifiers) {
  # extract gene names from TCGA RNASeq identifiers, use the entrez ID 
  # provided if no gene name is given. 
  # Ex: COL11A1|1301 -> COL11A1, ?|100133144 -> 100133144
  # 
  # Args:
  #   identifiers: identifiers provided by TCGA RNASeq matrix
  #
  # Returns:
  #   The gene names or the entriez IDs if gene names is `?`
  sapply(identifiers, function(x){
    tokens <- strsplit(x, split="\\|")[[1]]
    if (tokens[1] == "?"){
      return (tokens[2])
    } else {
      return (tokens[1])
    }
  })
}
