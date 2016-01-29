GSE51372_readCounts <- read.delim("~/Downloads/GSE51372/GSE51372_readCounts.txt", header=FALSE)
View(GSE51372_readCounts)

GSE51372 <- GSE51372_readCounts[,-(c(1,2,3,5,6))]
transformRead <- function(test){
  row.namesGSE5 <- as.character(unlist(test[,1]))[-1]
  col.namesGSE5 <- as.character(unlist(test[1,]))[-1]
  test <- test[,-1] 
  test <- test[-1,]
  new <- apply(as.matrix.noquote(test),2,as.numeric)
  rownames(new) <- row.namesGSE5
  colnames(new) <- col.namesGSE5
  #load charactersitics to look at sample data features
  return(new)
}