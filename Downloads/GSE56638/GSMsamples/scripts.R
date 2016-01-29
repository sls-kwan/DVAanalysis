setwd("~/IT/Downloads/GSE56638/GSMsamples/")
temp <- list.files(pattern = "*.txt")

myfiles = lapply(temp, function(x,y) read.delim(x, header = F))

test <- do.call(cbind, myfiles)

test <- test[,-seq(3, ncol(test), 2)]

row.namesGSE5 <- test$V1
test$V1 <- NULL
rawCountGSE56638 <- apply(as.matrix.noquote(test),2,as.numeric)
rownames(rawCountGSE56638) <- row.namesGSE5

#Getting the column names from the samples names
#We split up the filenames based on _ to get the raw GSM
tempUscore <- strsplit(temp, "_")
GSMGSE56638 <- lapply(tempUscore, function(x){return(x[1])})

GSMGSE56638 <- as.character(unlist(GSMGSE56638))

##Using the raw GSM when to use the Web Scraping tool GEOmetascrape.py
##To get the selected fields of the characteristics from each entry
## Needs two arguments Input file and Output File
tissues <- read.table("~/IT/Downloads/GSE56638/GSMsamples/tissues.tsv", quote="\"", comment.char="")
age <- read.csv("~/IT/Downloads/GSE56638/GSMsamples/age.tsv", header=FALSE, sep=";")
age <- lapply(age, function(x){gsub("developmental stage: ","", x)})$V1
age <- as.character(unlist(lapply(age, function(x){gsub("embryonic day ", "E", x)})))
colnames(rawCountGSE56638) <- mapply(function(x,y,z){return(paste(y,".", z, ".",x,sep=""))}, 1:length(tissues$V2), tissues$V2, age)

#19 of each, cardio and neuronal in that order
save(rawCountGSE56638, file="rawCountGSE56638.RData")
