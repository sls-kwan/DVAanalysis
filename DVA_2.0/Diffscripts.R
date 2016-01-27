source("~/IT/DVA_2.0/DiffAnalysis.R")
require(biomaRt)
DVAanalysis <- function(d, ngroup1, ngroup2, threshold, name){
  ##Need to include the pvalues as well
  dv <- DVAnalysis(d, ngroup1, ngroup2, threshold)
  de <- DEAnalysis(d, ngroup1, ngroup2, threshold)
  nonsignDV <- setdiff(rownames(d), dv)
  
  deNames <- list(rownames(de$UpDE), rownames(de$DownDE))
  nonsignDE <- setdiff(rownames(d), unlist(deNames))
  names(deNames) <- c("UpDE", "DownDE")
  both <- DEnDVAnalysis(dv, de)
  bothOpp <- DEUDVUAnalysis(dv, de)
  write(unlist(names(dv$UpDV)), file=paste("Results/Updv", name, ".txt", sep=""))
  write(unlist(names(dv$DownDV)), file=paste("Results/Downdv", name, ".txt", sep=""))
  
  write(rownames(de$UpDE), file=paste("Results/UpDE", name, ".txt", sep=""))
  write(rownames(de$DownDE), file=paste("Results/DownDE", name, ".txt", sep=""))
  
  write(unlist(both$Up), file=paste("Results/bothUp", name, ".txt",  sep=""))
  write(unlist(both$Down), file=paste("Results/bothDown", name, ".txt", sep=""))
  
  write(unlist(bothOpp$DownDVDEUp), file=paste("Results/DownDVDEUp", name, ".txt", sep=""))
  write(unlist(bothOpp$UpDVDEDown), file=paste("Results/UpDVDEDown", name, ".txt", sep=""))
  tableValues <- cbind(length(both$Down),length(intersect(deNames[[2]], nonsignDV)),length(intersect(deNames[[2]], dv$UpDV)),
                       length(intersect(dv$DownDV, nonsignDE)),length(intersect(nonsignDE, nonsignDV)),length(intersect(dv$UpDV, nonsignDE)),
                       length(intersect(dv$DownDV, deNames[[1]])),length(intersect(deNames[[1]], nonsignDV)),length(both$Up)
                       )
  tableValues <- matrix(tableValues, nrow = 3, ncol =3, byrow = T)
  write.table(tableValues, file=paste("Results/tableValues_", name, ".txt", sep=""), row.names = F, col.names = F)
  return(list(dv, deNames, both, bothOpp))
}

load("cancVnoncancFinal.RData")
d.control <- cancVnoncanc.final[,grep("H358", colnames(cancVnoncanc.final))]
d.cancer <- cancVnoncanc.final[,grep("LC-PT-45_", colnames(cancVnoncanc.final))]
d <- cbind(d.cancer, d.control)
rownames(d) <- gsub(".\\d+$", "", rownames(d))
resultsGSE69405 <- DVAanalysis(d, ncol(d.cancer), ncol(d.control), 0.01, "GSE69405")
save(resultsGSE69405, file="~/IT/DVA_2.0/resultsGSE69405.RData")

#load("nsmbDraft.RData")
#http://www.nature.com/nsmb/journal/v20/n9/full/nsmb.2660.html#supplementary-information

load("~/IT/Downloads/GSE56638/GSMsamples/rawCountGSE56638.RData")

##Please check column names to check for the see reasoning
##The first 19 samples are E14.5 and remaining 12 are adult cells
GSE56638.cardio <- rawCountGSE56638[,c(1:19, 39:50)]

##These are comparison of the cortex of brain and heart
##First 19 are heart and 19 are neuronal
GSE56638.heartNeuronal <- rawCountGSE56638[,-c(39:50)]

resultsGSE56638 <- DVAanalysis(GSE56638.heartNeuronal, 19, 19, 0.1, "GSE56638")
resultsGSE56638.heart <- DVAanalysis(GSE56638.cardio, 19, 12, 0.1, "GSE56638_heart")
save(resultsGSE56638, file="~/IT/Downloads/GSE56638/GSMsamples/resultsGSE56638.RData")
save(GSE56638.heartNeuronal,file="~/IT/Downloads/GSE56638/GSMsamples/GSE56638heartNeuronal.RData")
save(GSE56638.cardio, file="~/IT/Downloads/GSE56638/GSMsamples/GSE56638cardio.RData")

load("~/IT/Downloads/E-MTAB-2805.processed/rawCountCell.RData")
##This dataset has three groups S,G1, G2M
rawCountCell <- rawCountCell[-(38385:nrow(rawCountCell)),]
cellS <- rawCountCell[,grep("S_cell", colnames(rawCountCell))]
cellG1 <- rawCountCell[,grep("G1_cell", colnames(rawCountCell))]
cellG2M <- rawCountCell[,grep("G2M_cell", colnames(rawCountCell))]
cell <- cbind(cellG1, cellS)
save(cell, file="~/IT/Downloads/E-MTAB-2805.processed/CellG1S.RData")
# cellG1S <- cbind(cellG1, cellS)
# cellG1S <- parser(cellG1S)
#cell <- QCparser(cell)
resultsEMATB2805 <- DVAanalysis(cell, ncol(cellG1), ncol(cellS), 0.01, "E-MTAB-2805")
save(resultsEMATB2805, file="~/IT/Downloads/E-MTAB-2805.processed/resultsEMATB2805G1S.RData")
cellUp <- as.numeric(unlist(lapply(resultsEMATB2805G1S[[3]]$Up, function(x){which(x == rownames(cell))})))
load("~/IT/Downloads/GSE67835/GSE67835.RData")
GSE67836.AB_S5 <- GSE67835[,grep("AB_S5", colnames(GSE67835))]
GSE67835.FB_S2 <- GSE67835[,grep("FB_S2", colnames(GSE67835))]

GSE67836.final <- cbind(GSE67836.AB_S5, GSE67835.FB_S2)

rownames(GSE67836.final) <- gsub("\\s", "", unlist(rownames(GSE67836.final)))
GSE67836.final <- QCparser(GSE67836.final)
save(GSE67836.final, file="~/IT/Downloads/GSE67835final.RData")
resultsGSE67836 <- DVAanalysis(GSE67836.final, ncol(GSE67836.AB_S5), ncol(GSE67835.FB_S2), 0.01, "GSE67835")
save(resultsGSE67836, file="~/IT/Downloads/GSE67835/resultsGSE67836.RData")
# load("~/IT/Downloads/GSE45719/GSE45719.RData")
# 
# GSE45719_16cell <- GSE45719[,1:50]
# GSE45719_8cell <- GSE45719[,51:64]
# 
# d <- cbind(GSE45719_16cell, GSE45719_8cell)
# 
# dv <- DVAnalysis(d, 50, 14)
# de <- DEAnalysis(d, 50, 14)
# 
# both <- DEnDVAnalysis(dv, de)
# mostdvGSE45179 <- gsub("\\s", "", unlist(dv))
# write(mostdvGSE45179, file="GSE45179.txt")


load("~/IT/Downloads/GSE51372/GSE51372.RData")
##Note that more multiple rownames limma will produce GeneIds in ID column rather
##than in rownames and this has multiple gene names that are the same

##This can be remedied using the maximum median approach, selecting
##the gene isoform with highest gene name
d.cancer <- GSE51372[,grep("TuMP3.*", colnames(GSE51372))]
d.cancerControl <- GSE51372[,grep("nb508.*", colnames(GSE51372))]
d <- cbind(d.cancer, d.cancerControl)
d <- QCparser(d)
##Remove novel genes
d <- d[-(which(rownames(d) == "")),]
duplicatedRowNames <- rownames(d)[which(duplicated(rownames(d)))]
duplicatedRowIndex <- which(duplicated(rownames(d)))

#Function that generates the indexs of the non max median 
#duplicated rows, so we can exclude them
notMaxMedian <- lapply(unique(duplicatedRowNames), function(x){
  duplRows <- d[which(x == rownames(d)),]
  duplIndex <- which(x == rownames(d))
  medianValues <- apply(duplRows, 1, function(x){
    return(median(x))
  })
  return(duplIndex[-which.max(medianValues)])
})

d <- d[-(unlist(notMaxMedian)),]
resultsGSE51372 <- DVAanalysis(d, 16, 9, 0.1, "GSE51372")
save(resultsGSE51372, file="~/IT/Downloads/GSE51372/resultsGSE51372.RData")
# load("~/IT/Downloads/GSE75104/GSE75104.RData")
# resultsGSE75104 <- DVAanalysis(QCparser(GSE75104), 15, 15, 0.1, "GSE75104")
# 
# load("~/IT/Downloads/GSE75108/GSE75108.RData")
# load("~/IT/Downloads/GSE75107/GSE75107.RData")
# GSE7510x <- cbind(GSE75108, GSE75107)
# resultsGSE7510x <- DVAanalysis(QCparser(GSE7510x), ncol(GSE75108), ncol(GSE75107), 0.01, "GSE7510x")
