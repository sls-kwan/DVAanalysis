##Remove zero 0.95
QCparser <- function(d){
  PercentZero<-apply(d,1, function(x){
    length(which(x==0))/length(x)
  });
  return(d[which(PercentZero <= 0.95),]);
}
##Creates plot for each RNA-seq sample
displayExpression <- function(d, name){
  setwd("~/IT/DVA_2.0/QC/")
  system(paste("mkdir", name))
  for (i in 1:ncol(d)){
      if(i %% 9 == 0){
      plot(density(log2(d[,i]+1))); abline(v=2)
      dev.off()
      if((i/9)+1 < 10){
        png(filename=paste(name, "/",name, "_", "0",(i/9)+1, ".png", sep=""))
      } else{
        png(filename=paste(name, "/", name, "_", (i/9)+1, ".png", sep=""))
      }
      par(mfrow=c(3,3))
      next()
      } else if(i == 1){
        png(filename=paste(name, "/", name, "_", "0",i, ".png", sep=""))
        par(mfrow=c(3,3))
      } else{}
      plot(density(log2(d[,i]+1))); abline(v=2)
  }
  dev.off()
}

LookedAt <- read.csv("~/IT/DVA_2.0/LookedAt.csv", header=FALSE, comment.char ="#")
images <- apply(LookedAt, 1, function(x){
  data <- load(as.character(x[2]))
  storeData <- get(data)
  displayExpression(storeData, as.character(x[1]))
  system(paste("./convertPngPdf.sh",x[1]))
})
nrowD <- 1:nrow(cell)
splitRefer <- split(cellUpValues, ceiling(seq_along(cellUpValues)/10))
splitRefer <- split(cell, ceiling(seq_along(cell)/10))
name <- c("GSE56638")

##Use all of data, but subset it into selecting the ones that are cellUp or cellDown and split on those
mapply(function(x, y){
  if((y/10)+1 < 10){
    png(filename=paste("bothGO", "/",name, "/", name, "_", "0",y, ".png", sep=""))
  } else{
    png(filename=paste("bothGO", "/", name, "/", name,"_", y, ".png", sep=""))
  }
visDVList(cell,x, n=c(ncol(cellG1), ncol(cellS)), rank=FALSE,sorted=TRUE,label=c("Heart","Neuronal"),nCols=5, nRow=-1, title=paste(rownames(cell)[x]))
dev.off()
}, splitRefer, as.numeric(names(splitRefer)))

cellDown <- match(resultsGSE51372[[3]]$Down, rownames(d))
name <- "GSE51372"
name2 <- "GSE51372Down"
mapply(function(x, y){
  if((y/10)+1 < 10){
    png(filename=paste("bothGO/", name, "/", name2, "_", "0",y, ".png", sep=""))
  } else{
    png(filename=paste("bothGO/", name, "/", name2, "_", y, ".png", sep=""))
  }
  visDVList(d,x, n=c(16,9), rank=FALSE,sorted=TRUE,label=c("Cancer","CTCs"),nCols=5, nRow=-1, title=paste(rownames(d)[x]))
  dev.off()
# }, splitRefer, as.numeric(names(splitRefer)))
}, splitReferD, as.numeric(names(splitReferD)))
splitRefer <- split(cellUp, ceiling(seq_along(cellUp)/10))
splitReferD <- split(cellDown, ceiling(seq_along(cellDown)/10))
