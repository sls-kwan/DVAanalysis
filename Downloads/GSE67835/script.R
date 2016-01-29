setwd("~/IT/Downloads/GSE67835/")
temp <- list.files(pattern = "*.csv")

myfiles = lapply(temp, function(x,y) read.delim(x, header = F))

test <- do.call(cbind, myfiles)

test <- test[,c(1,seq(2,ncol(test),2))]

test <- test[1:(nrow(test)-3),]

row.namesGSE5 <- as.character(unlist(test$V1))
test$V1 <- NULL
GSE67835 <- apply(as.matrix.noquote(test),2,as.numeric)
rownames(GSE67835) <- row.namesGSE5
load("~/IT/Downloads/GSE67835/GSMmetadata.RData")
lengthGSM <- c(1:length(GSMmetadata))
colNames <- lapply(lengthGSM, function(x){
  paste(GSMmetadata[x], ".", x, sep="")
})
colnames(GSE67835) <- colNames
#19 of each, cardio and neuronal in that order


#GSE67835.test <- GSE67835[,order(GSMmetadata)]
#load("~/Downloads/GSE67835/GSMmetadata.RData")
#colnames(GSE67835.test) <- as.character(unlist(sort(GSMmetadata)))
save(GSE67835, file="GSE67835.RData")
