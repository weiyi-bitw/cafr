IsCancerTCGA <- function(x) {
  # check if a TCGA sample ID is a cancer sample or not
  #
  # Args:
  #   x: TCGA sample ID (TCGA-XX-XXXX-XX)
  #
  # Returns:
  #   TRUE if this sample is a cancer sample
  #
  substr(x, 14, 14) == "0"
}
