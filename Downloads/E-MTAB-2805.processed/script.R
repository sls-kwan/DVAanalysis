setwd("~/Downloads/CellCycle/")
temp <- list.files(pattern = "*.txt")

myfiles = lapply(temp, function(x,y) read.delim(x, header = F))

test <- do.call(cbind, myfiles)

test <- test[,-which(colnames(test) == c("V1", "V2", "V3", "V4"))[-1]]

row.namesGSE5 <- as.character(unlist(test$V1))
test$V1 <- NULL
col.namesGSE5 <- as.character(unlist(test[1,]))
test <- test[-1,]
rawCountCell <- apply(as.matrix.noquote(test),2,as.numeric)
rownames(rawCountCell) <- row.namesGSE5[-1]
colnames(rawCountCell) <- col.namesGSE5
dontNeed <- rownames(tail(rawCountCell))
rawCountCell <- rawCountCell[-which(rownames(rawCountCell) == dontNeed),]
rawCountCell <- rawCountCell[-grep("ERCC", rownames(rawCountCell)),]
#19 of each, cardio and neuronal in that order
save(rawCountCell, file="rawCountCell.RData")
