stringsAsFactors=F

source("~/genequery/R/genequery_functions.R")

args <- commandArgs(TRUE)
csv  <- args[1] 

library(BiocParallel) 
library (WGCNA)

t0 <- proc.time()
## no parallel execution in this one. 

parsens_single_wgcna(csv)

proc.time()-t0
