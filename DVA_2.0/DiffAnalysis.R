setwd("~/IT/DVA_2.0/")
source("DVA_2_0Re.R")

##Input, d is dataframe, group1 and group2 are numbers in each respective group
##Cutoff is the minimum cutoff
DVAnalysis <- function(d, group1, group2 , cutoff){
  ##Creates DV analysis Framework
  result<-dva(d,group1,group2,test="F",outlier.rm=T,r=1.5,piter=5000,thres=0,verbose=F);
  
  ##Adjusts the P value for all values due multiple testing
  adjustPvalue <- p.adjust(result$p.dv, method="fdr", n = nrow(d))
  
  #Selects elements that satisfies the cutoffs and also above zero
  pThreshold <- which(adjustPvalue <= cutoff & adjustPvalue > 0.000000e+00)
  
  #Get the coressponding names and order them and return also their pvalues
  names.pThreshold <- rownames(d)[pThreshold]
  opThresol <- order(adjustPvalue[pThreshold])
  indexDVGenes <- opThresol
  pValues <- adjustPvalue[indexDVGenes]
  
  #hist(result$p.dv[which(adjustPvalue <= 0.01 & adjustPvalue > 0.005)], xlab="p",main="distribution of p value for DV");
  #Gets the Pvalues and gene names and then uses two different groups
  #DV using variance and then using the difference betwen variances > 0 
  #To determine if something is up or down
  DVGenes<- names.pThreshold[indexDVGenes]
  control <- d[indexDVGenes, c(0:group1)]
  disease <- d[indexDVGenes, -c(0:group1)]
  
  controlVar <- apply(control, 1, function(x){return(var(x, na.rm=TRUE))})
  diseaseVar <- apply(disease, 1, function(x){return(var(x, na.rm=TRUE))})
  
  upDVGenes <- DVGenes[which(diseaseVar-controlVar >0)]
  upPValues <- pValues[which(diseaseVar-controlVar >0)]
  downDVGenes <- DVGenes[which(diseaseVar-controlVar <0)]
  DownpValues <- pValues[which(diseaseVar-controlVar <0)]
  names(upPValues) <- upDVGenes
  names(DownpValues) <- downDVGenes
  allDVGenes <- list(upPValues, DownpValues)
  names(allDVGenes) <- c("UpDV", "DownDV")
  return(allDVGenes)
}

##Input, d is dataframe, group1 and group2 are numbers in each respective group
##Threshold is the minimum cutoff
DEAnalysis <- function(d, group1 , group2, threshold){
  ##Differential Expression analysis using the limma and edge R package
  require(limma)
  require(edgeR)
  
  #Creating the different cell types, creating model matrix and doing
  #voom normalization on data
  cellType <- factor(c(rep("cancer", group1), rep("control", group2)))
  
  design <- model.matrix(~cellType)
  
  dge <- DGEList(counts=d)
  
  dge <- calcNormFactors(dge)
  
  v <- voom(dge,design,plot=F)
  
  fit <- lmFit(v,design)
  
  fit <- eBayes(fit)
  
  #Creating a topTable of the top genes that have a lower value than the threshold 
  #High number because limma package default value is 10
  DEgenes <- topTable(fit,coef=ncol(design), number = 100000000, p.value = threshold)
  
  UpRegulated <- DEgenes[which(DEgenes$logFC > 0),]
  
  DownRegulated <- DEgenes[which(DEgenes$logFC < 0),]
  allDEgenes <- list(UpRegulated, DownRegulated)
  names(allDEgenes) <- c("UpDE", "DownDE")
  return(allDEgenes)
}

##Function that finds Both Up and Down in DE and DV analysis
## Uses a list of DVgenes and DEgenes
DEnDVAnalysis <- function(DVgenes, DEgenes){
  UpDVnDE <- intersect(names(DVgenes$UpDV), rownames(DEgenes$UpDE))
  DownDVnDE <- intersect(names(DVgenes$DownDV), rownames(DEgenes$DownDE))
  UpnDown <- list(UpDVnDE, DownDVnDE)
  names(UpnDown) <- c("Up", "Down")
  return(UpnDown)
}

##Function that finds Up and Down (Overlapping) DV and DE analysis in gene list
##Uses a list of DVgenes and DEgenes
DEUDVUAnalysis <- function(DVgenes, DEgenes){
  DownDVnDE <- intersect(names(DVgenes$UpDV), rownames(DEgenes$DownDE))
  UpDVnDE <- intersect(names(DVgenes$DownDV), rownames(DEgenes$UpDE))
  DVnDE <- list(DownDVnDE, UpDVnDE)
  names(DVnDE) <- c("DownDVDEUp", "UpDVDEDown")
  return(DVnDE)
}



