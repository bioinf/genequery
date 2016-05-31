wgcna_cluster = function(csvfile,ngenes=6000) {
     
  cat(sprintf("Processing file %s, top %d most expressed genes will be considered\n",csvfile,ngenes))
  
  ## this assumes the pre-processed format that's output by preprocess_matrix.pl 
  tt <- read.csv(csvfile,row.names=1)
  ttf <- tt[grep("NONE",tt$Symbol,invert=T),]
  ttf <- ttf[grep("NONE",ttf$Entrez_ID,invert=T),]
  ttc <- collapseRows(ttf[3:ncol(ttf)],ttf$Entrez_ID,rownames(ttf))
  
  cat(sprintf("Initial expression matrix dimensions:\n"))
  cat(sprintf("%d\n",dim(tt)))
  cat(sprintf("Defined probe-only expression matrix dimensions:\n"))
  cat(sprintf("%d\n",dim(ttf)))
  cat(sprintf("Collapsed expression matrix dimensions:\n"))
  cat(sprintf("%d\n",dim(ttc$datETcollapsed)))
  
  exp <- as.data.frame(ttc$datETcollapsed)
  exp$max2 <- apply(exp,1,FUN=max2)
  exp <- exp[order(exp$max2,decreasing=T),]
  exp <- exp[1:ngenes,]
  exp$max2 <- NULL 
  cat(sprintf("Truncated collapsed expression matrix dimensions:\n"))
  cat(sprintf("%d\n",dim(exp)))
  
  texp <- t(exp)
  pwrs <- c(5:25)
  cat(sprintf("======================================================================\n"))
  sft1 <- pickSoftThreshold(texp,powerVector=pwrs,verbose=3,networkType="signed hybrid",RsquaredCut = 0.7)
  sft2 <- pickSoftThreshold(texp,powerVector=pwrs,verbose=3,networkType="signed hybrid",RsquaredCut = 0.8)
  sft3 <- pickSoftThreshold(texp,powerVector=pwrs,verbose=3,networkType="signed hybrid",RsquaredCut = 0.9)
  cat(sprintf("======================================================================\n"))
  TAG <- gsub("_preprocessed.csv","",csvfile)
  
  p1 <- sft1$powerEstimate
  p2 <- sft2$powerEstimate
  p3 <- sft3$powerEstimate
  
  if (!is.na(p3)) {
    ncut <- if(p3 <= 10) 0.25 else 0.15
    cat(sprintf("Performing clustering for mergeCutHeight %s and power %d\n",ncut,p3))
    net <- blockwiseModules(texp, power = p3, minModuleSize = 50, reassignThreshold = 0, mergeCutHeight = ncut,
                            maxBlockSize=6000,deepSplit=4,numericLabels=T,pamRespectsDendro=F,networkType = "signed hybrid",verbose = 0)
    print(table(net$colors))
    write.table(cbind(rownames(exp),net$colors),file=paste(TAG,"_hybrid_modules.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
    write.table(t(net$MEs),file=paste(TAG,"_eigengenes.txt",sep=""),quote=F,col.names=colnames(exp),row.names=T,sep="\t")

  } else if (!is.na(p2)) {
    ncut <- if(p2 <= 10) 0.25 else 0.15
    cat(sprintf("Performing clustering for mergeCutHeight %s and power %d\n",ncut,p2))
    net <- blockwiseModules(texp, power = p2, minModuleSize = 50, reassignThreshold = 0, mergeCutHeight = ncut,
                            maxBlockSize=6000,deepSplit=4,numericLabels=T,pamRespectsDendro=F,networkType = "signed hybrid",verbose = 0)
    print(table(net$colors))
    write.table(cbind(rownames(exp),net$colors),file=paste(TAG,"_hybrid_modules.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
    write.table(t(net$MEs),file=paste(TAG,"_eigengenes.txt",sep=""),quote=F,col.names=colnames(exp),row.names=T,sep="\t")
    
  } else if (!is.na(p1)) {
    ncut <- if(p1 <= 10) 0.25 else 0.15
    cat(sprintf("Performing clustering for mergeCutHeight %s and power %d\n",ncut,p1))
    net <- blockwiseModules(texp, power = p1, minModuleSize = 50, reassignThreshold = 0, mergeCutHeight = ncut,
                            maxBlockSize=6000,deepSplit=4,numericLabels=T,pamRespectsDendro=F,networkType = "signed hybrid",verbose = 0)
    print(table(net$colors))
    write.table(cbind(rownames(exp),net$colors),file=paste(TAG,"_hybrid_modules.txt",sep=""),quote=F,col.names=F,row.names=F,sep="\t")
    write.table(t(net$MEs),file=paste(TAG,"_eigengenes.txt",sep=""),quote=F,col.names=colnames(exp),row.names=T,sep="\t")
    
  } else {
    cat(sprintf("No appropriate clustering conditions were found!\n"))
  }
  
  cat(sprintf("CLUSTERING IS COMPLETE. See files *modules.txt for the list of modules and  *eigengenes.txt for averaged module expression.\n"))
}

max2 = function(array) { 
	n = length(array)
	sort(array,partial=n-1)[n-1]
}

make_svg_heatmaps = function (filename) {
  ## With non-renamed/non-normalized eigengene tables.
  ##library(ggplot2)
  ##library(reshape)
  ##rm(list=ls(all=TRUE)) ## crucial here!
  tag <- gsub("_eigengenes.txt","",filename)
  tt <- read.table(filename,header=T,row.names=1)
  row.names(tt) <- as.numeric(gsub("ME","",row.names(tt)))
  tts <- tt[order(as.numeric(row.names(tt))),]
  tts2 <- tts[,order(colnames(tts))]
  ttn <- as.data.frame(t(apply(tts2,1,FUN=exvector)))
  ttn$Module <- row.names(ttn)
  
  nsmp <- ncol(tt)
  nmod <- nrow(tt)-1 ## number of non-zero modules
  tt.m <- melt(ttn)
  
  names(tt.m) <- c("Module","Sample","value")
  ## the following are very important
  tt.m$Module <- as.factor(tt.m$Module)
  sorted_labels <- paste(sort(as.integer(levels(tt.m$Module))))
  tt.m$Module <- factor(tt.m$Module, levels = sorted_labels)
  tt.m$Sample <- with(tt.m, factor(Sample,levels = rev(sort(unique(Sample)))))
  
  ## let's evaluate SVG dimensions. We want constant width of 10 in. 
  title_length <- max(nchar(names(tt)))
  svg_height  <- round((10/27.1)*nsmp,digits=3)
  svg_width <- round((10/27.1)*(0.25*title_length+nmod+1),digits=3)
  cat (sprintf("Longest sample title is %d characters, there are %d samples and %d modules, plot width is %f, height is %f\n",title_length,nsmp,nmod,svg_width,svg_height))
  
  base_size <- 9
  for (i in 0:nmod) {
    ## special type of coloring for easier heatmap reading. Can be adjusted. 
    p <- ggplot(tt.m, aes(Module,Sample)) + geom_tile(aes(fill = value), colour = "black",size=0.5) + scale_fill_gradientn(colours=c("blue","blue","white","red","red"),space="Lab")
    p2 <- p + theme_grey(base_size=base_size)+labs(x="",y="")+scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand =c(0,0))+theme(legend.position ="none",axis.text.y=element_text(size=base_size*1.2,colour="black"),axis.text.x=element_text(size=base_size*1.2,colour="black"))
    xmin <- i+0.5
    xmax <- i+1.5
    ymin <- 0.5
    ymax <- nsmp+0.5
    rect <- data.frame(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    p3 <- p2 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),color="black",fill="grey",alpha=0, size=3, inherit.aes = F)
    svgname <- paste(tag,"_module_",i,".svg",sep="")
    ##cat(sprintf("Processing module number %d, saving file %s\n",i,svgname))
    ggsave(file=svgname, plot=p3, width=svg_width, height=svg_height,limitsize=F)
  }
}

exvector = function (v) {
  min=min(v)
  max=max(v)
  for (i in 1:length(v)) {
    if (v[i] <= 0) {
      v[i]=round(v[i]/abs(min),digits=3)
    } else { 
      v[i]=round(v[i]/max,digits=3)
    }  
  }
  return(v)
}
