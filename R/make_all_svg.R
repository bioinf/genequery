source("~/genequery/R/genequery_functions.R")

## run this on 770, there's some version discrepancy which makes defective SVGs on 769
## run this in the directory with all the eigengene text files 

library(ggplot2)
library(reshape)
library(svglite)
stringsAsFactors=F

lf=list.files(pattern="eigengenes.txt")

cat(sprintf("Found %d files with the appropriate extension. Begin processing...\n",length(lf)))

for (i in 1:length(lf)) {
  cat(sprintf("Processing file number %d.\n",i))
  make_svg_heatmaps(lf[i])
}

cat(sprintf("FINISHED generating the plots. Bye now!\n"))
