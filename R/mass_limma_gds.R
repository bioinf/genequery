mass_limma_gds = function (tag) {
  library(limma)
  library(WGCNA)
  library(data.table)
  
  exp_file   <- paste(tag,"exp",sep=".")
  cond_file  <- paste(tag,"cond",sep=".")
  pairs_file <- paste(tag,"contrast",sep=".")
  
  ## there's a stupid probe name with a #, it gets disabled with comment.char=""
  tt         <- read.table(exp_file,comment.char="",header=T,row.names=1,sep="\t",fill=T,check.names=F) 
  tt         <- tt[complete.cases(tt),]
  
  tt2        <- collapseRows(tt[2:ncol(tt)],tt$IDENTIFIER,rownames(tt))
  tt3        <- as.data.frame(tt2$datETcollapsed)
  
  if (mean(colMeans(tt3))>12) {
    tt3 <- log2(tt3) 
  }
  tt3        <- tt3[complete.cases(tt3),]
  
  tt3$max2   <- apply(tt3,1,FUN=max2)
  ngenes     <- min(nrow(tt3),6000)
  tt4        <- tt3[order(tt3$max2,decreasing=T),]
  tt4        <- tt4[1:ngenes,]
  tt4$max2   <- NULL 
  
  exp_table  <- tt4
  cond_table <- read.table(cond_file,header=T,row.names=1)
  
  DE       <- list()
  pairs    <- fread(pairs_file) ## data.table object with needed column, cond1, and cond2
  N        <- nrow(pairs)

  for (i in 1:N) {
    tryCatch({
      column <- pairs$Column[i]
      cond1  <- pairs$cond1[i]
      cond2  <- pairs$cond2[i]
      cat(sprintf("GDS-LIMMA: Index: %d, cond1: %s, cond2 %s\n",i,cond1,cond2))
      DE[[i]] <- limma_gds(exp_table,cond_table,column,cond1,cond2)
    }, error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
  }
  return(DE)
}
