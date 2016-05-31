source("~/genequery/functions.R")

## run this in the directory that contains the *eigengenes.txt files obtained after WGCNA clustering. 
library(ggplot2)
library(reshape)
stringsAsFactors=F

lf=list.files(pattern="eigengenes.txt")

cat(sprintf("Found %d files with the appropriate extension. Begin processing...\n",length(lf)))

for (i in 1:length(lf)) {
  cat(sprintf("Processing file number %d.\n",i))
  make_svg_heatmaps(lf[i])
}

cat(sprintf("FINISHED generating the plots. Bye now!\n"))
