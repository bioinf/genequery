source("~/genequery/functions.R")

args <- commandArgs(TRUE)
csv  <- args[1] 

if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      A script for WGCNA production clustering.
 
      Example:
      ./wgcna_cluster.R preprocessed.csv \n\n")
 
  q(save="no")
}

library (WGCNA)
stringsAsFactors=F

wgcna_cluster(csv)

