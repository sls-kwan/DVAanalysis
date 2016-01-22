##################################################################################
#  DVA: Differential Variability Analysis of Gene Expression in R                #
#  Copyright (C) 2013 Joshua W.K. Ho (j.ho @AT@ victorchang .DOT. edu .DOT. au)  #
#                                                                                #
#    This program is free software; you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by        #
#    the Free Software Foundation; either version 2 of the License, or           #
#    (at your option) any later version.                                         #
#                                                                                #
#    This program is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of              #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
#    GNU General Public License for more details.                                #
#                                                                                #
#    You should have received a copy of the GNU General Public License along     #
#    with this program; if not, write to the Free Software Foundation, Inc.,     #
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                 #
#                                                                                #
#                                                                                #
##################################################################################


###
# For Visualization
###

###
# A general function for plotting n groups of expression value of a single gene across k samples
# geneExp: a vector of gene expression value across a number of arrays
# n: an array that specify the number of samples in each group for example: n=c(10,20,30) for three groups with 10,20 and 30 arrays each
# rank: a boolean to specify if the value should be ranked, default is false
# sorted: a boolean to specify if expression value in each group should be sorted
# order: a vector to specify the ordering of arrays (when sorted is FALSE)
# label: a vector of label of each group
# title: the title of the graph
###
visDV<-function(geneExp, n=c(), rank=FALSE, sorted=TRUE, order=1:length(geneExp), label=n, title="Gene"){
   if(rank){
      geneExp<-rank(geneExp);
   }
   numArray<-length(geneExp);
   numGroups<-length(n);
   if(sum(n) != numArray){
      print("ERROR: the total number of array in n does not match the number of array in geneExp");
      return("Error");
   }
   if(length(label) != numGroups){
      label <-rep("", numGroups);
   }
   pt = 1;  #A pointer to the current index
   groupAll<-c();
   if(sorted){
      for(i in n){
         endPt <- pt + i - 1;
         groupi <- sort(geneExp[pt:endPt]);
         groupAll <- c(groupAll, groupi);
         pt <- endPt + 1;
      }
   }
   else{
      groupAll<-geneExp[order];
   }
   plot(groupAll, main=title, xlab="array (sorted)", ylab="expression");
   pt <- 1;
   midY<-max(groupAll);
   j<-1;
   for(i in n){
      endPt <- pt + i - 1;
      midPt <- (pt + endPt)/2;
      text(midPt, midY, label[j]);
      if(endPt != numArray){
         abline(v=(endPt+0.5));   #draw the line separating the groups
      }
      pt <- endPt + 1;
      j <- j+1;
   }
   
   return(groupAll);
}

###
# A general function for plotting n groups of expression value of a list of genes across k samples
# data: a data frame containing the raw data
# geneList: an array of indices of genes to be analyzed
# n: an array that specify the number of arrays in each group for example: n=c(10,20,30) for three groups with 10,20 and 30 arrays each
# rank: a boolean to specify if the value should be ranked, default is false
# sorted: a boolean to specify if expression value in each group should be sorted
# label: a vector of label of each group
# nCols: the number of columns per row
# nRow: the number of row 
# title: the title of each graph (for each gene)
###
visDVList<-function(data,geneList,n=c(),rank=FALSE,sorted=TRUE,label=n,nCols=5, nRow=-1, title=rep("title",length(geneList))){
   nGene<-length(geneList);
   nArray<-length(data[1,]);
   nCol<-nCols;
   if(nRow == -1){
      nRow<-ceiling(nGene/nCol);
   }
   par(mfrow=c(nRow,nCol));
   for(i in 1:nGene){
      visDV(as.vector(data[geneList[i],]), n, rank,sorted,1:nArray,label,title[i]);
   }
}

###
# A function for plotting heatmap for differential expression analysis. The colour scheme is topo.colors(100).
# data: a data frame containing the raw data
# geneList: an array of indices of genes to be analyzed
# n: an array that specify the number of arrays in each group for example: n=c(10,20,30) for three groups with 10,20 and 30 arrays each
# title: the title of each graph (for each gene)
# clustering: indicate if clustering of row and column is to be performed
###
visHeatMap<-function(data, geneList, n=c(), title="", clustering=FALSE){
   numArray<-length(data[geneList[1],]);
   if(sum(n) != numArray){
      print("ERROR: the total number of array in n does not match the number of array in data");
      return("Error");
   }
   colArray<-c("navy", "green", "blue","yellow", "red","orange","grey","black","pink","navy");
   cols<-c();
   j<-1;
   for(i in n){
      cols<-c(cols, rep(colArray[j], i));
      j<-j+1;
   }
   if(clustering){
        heatmap(data[geneList,], scale="none",  col=topo.colors(100), ColSideColors=cols, main=title);
   }
   else{
        heatmap(data[geneList,], Rowv=NA, Colv=NA, scale="none",  col=topo.colors(100), ColSideColors=cols, main=title);
        # col = cm.colors(256)
   }
}


##### Examples #####
#visDV(sample, c(40,20,40), sorted=FALSE);
#visDVList(data,1:5, c(10,40,30,20), sorted=F, label=lab, nCol=3);
#visHeatMap(mat, c(1,2), c(10,10), title="hello")
##### END Examples #####

###
# For Power calculation
###

###
# A power function of the F-test
#  n1 = sample size for group 1
#  n2 = sample size for group 2
#  sig = significance level of the test
#  a = The true ratio between the two variances
###
power.F.test<-function(n1, n2, sig, a){
   pow<- -1;
   if(a>1){
      f.null <- qf(sig/2, n1-1, n2-1, 0, lower.tail=FALSE);
      pow<-pf(f.null/a, n1-1, n2-1, 0, lower.tail=FALSE);
   }
   else{
      f.null <- qf(sig/2, n1-1, n2-1, 0, lower.tail=TRUE);
      pow<-pf(f.null/a, n1-1, n2-1, 0, lower.tail=TRUE);
   }
   return(pow);
}

###
# For DV analysis
###

###
# Return a vector of values with outlier marked as a tag. Outlier removal is based on inter-quartile range criteria
#  data = a vector of expression value of one gene
#  r = multiplier of the inter-quartile range. Small r can identify more values as outliers
#  tag = a number that is placed into the vactor to replace an outlier   (default=-1)
#  verbose = a boolean to determine if it anything should be printed.
###
detectOutlier<-function(data,r=1.5,tag=NA,verbose=FALSE){
   size<-length(data);
   q<-quantile(data, var=TRUE);
   IQR<-q[4]-q[2];
   if(verbose){
      print(paste("IQR: ", IQR));
      print(paste("Q2: ", q[2]));
      print(paste("Q4: ", q[4]));
   }
   removeOutlier <- function(x, q, IQR,tag){
     if(x>q[4]+r*IQR || x<q[2]-r*IQR){
       x <- tag
     }
     return(x)
   }
   data <- sapply(data, removeOutlier,q, IQR,tag, simplify = T)
   return(data)
}

###
# Determine if a statement satisfy a criteria
#  series = a vector of numerical values
#  thres = a threshold 
#  test = test whether it is more or less than
###
pass<-function(series,thres,test="less"){
   if(test=="less"){
      length(which(series<thres))
   }
   else if(test=="more"){
      length(which(series>thres))
   }
   else if(test=="equal"){
      length(which(series==thres))
   }
   else{
      print("ERROR, pass() can only have parameter test being less or more");
      return("ERROR, pass() can only have parameter test being less or more");
   }
}


###
# Calculate robust variance estimate. It calculates Qn = [|xi-xj|,i<j](k). sd = 2.2219 MAD approximately
# This is only a straight-forward implementation, Rousseeuw and Croux (1993) has provided a smarter and faster implementation
#  x = a vector of numbers
###
Qn<-function(x){
  n<-length(x);
  all<-c();
  for(i in 1:n){
     for(j in  i:n){
        all<-c(all, abs(x[i]-x[j]));
     }
  }
  return(2.2219*sort(all)[choose(floor(n/2)+1,2)]);
}


###
# Perform a gene-by-gene test of differential variability and differential expression. It returns an (unadjusted) p-value
#   data = data matrix
#   df1 = number of samples in class 1
#   df2 = number of samples in class 2
#   test = the type of DV test, "F" for F-test by default
#   outlier.rm = removal of outlier arrays
#   r = the r parameter for outlier removal
#   piter = number of iteration in a permutation test
#   thres = lower threshold for expression
#   verbose = verbose mode?
###
dva<-function(data,df1,df2,test="F",outlier.rm=TRUE,r=1.5,piter=5000,thres=0,verbose=F){
   nGene<-length(data[,1]);
   nArrays<-length(data[1,]);
   #for error handling
   if(df1+df2!=nArrays){
      return("Error The number of column in the matrix is not df1+df2");
   }
   #go through each array element to perform a f-test
   p<-array(-1,nGene);
   actualDF1<-array(-1,nGene);
   actualDF2<-array(-1,nGene);
   statistic<-array(-1,nGene);
   t.value<-array(-1,nGene);
   p.dv<-array(-1,nGene);
   p.de<-array(-1,nGene);
   for(g in 1:nGene){
      d1<-df1;
      d2<-df2;
      if(verbose){
         print(paste(g," "));
      }
      exp<-as.vector(data[g,]);
      #See the amount of zeros that exceed the limit
      #Made redudant using QCparser
       if(pass(exp,thres,"more")<=1){  
          p.dv[g]<- 1;
          p.de[g]<- 1;
       }
       else{
         if(outlier.rm){
           expTemp<-detectOutlier(unlist(data[g,]),r);
           d1t<- d1-length(which(which(is.na(expTemp))<=d1));
           d2t<- d2-length(which(which(is.na(expTemp))>d1));
           d1<-d1t;
           d2<-d2t;
           exp<- expTemp[which(!is.na(expTemp))];
         }
         else{
           exp<-unlist(data[g,]);
         }
            actualDF1[g]<-d1;
            actualDF2[g]<-d2;
            range1<-1:d1;
            range2<-(d1+1):(d1+d2);
            exp1<-exp[range1];
            exp2<-exp[range2];
            if(test =="F"){
             var.exp1 <- var(exp1, na.rm = TRUE)
             var.exp2 <- var(exp2, na.rm = TRUE)
             ##This is to ensure when we do a two-tailed test that the pvalues are valid so that
             ##S = Var 0/0 does not occur (Set it to 0), we should set it to a p-value of 1
             if((var.exp1 == 0 | is.na(var.exp1))  && (var.exp2== 0 | is.na(var.exp2))){
                s <- 0
                p.dv[g] <- 1
             } else {
              s<- var(exp1, na.rm = TRUE)/var(exp2, na.rm = TRUE);
              p.dv[g]<-pf(s,d1-1,d2-1,lower.tail=FALSE);
             }
               #take the upper tail
             statistic[g]<-s;
            }
            else if(test =="MAD"){       #Robust estimator of based on median absolute deviation
             s<-rVar(exp1)/rVar(exp2)
             p.dv[g]<-pf(s,d1-1,d2-1,lower.tail=FALSE);  #take the upper tail
             statistic[g]<-s;
            }
            else if(test =="Qn"){     #Robust estimator of based on Qn             NOT STABLE YET!
               s<-(Qn(exp1)^2)/(Qn(exp2)^2);
               p.dv[g]<-pf(s,d1-1,d2-1,lower.tail=FALSE);  #take the upper tail
               statistic[g]<-s;
             }
            else if(test == "pF"){
                  s<- var(exp1)/var(exp2);
                  theta <- 0;
                  for(k in 1:piter){
                     randExp<-sample(exp);
                     rand1Temp<-randExp[range1];
                     rand2Temp<-randExp[range2];
                     varTemp1<-var(rand1Temp);
                     varTemp2<-var(rand2Temp);
                     if(varTemp1/varTemp2>s){
                        theta <-theta +1;
                     }
                  }
                  statistic[g]<-s;
                  p.dv[g]<-theta /piter;
             }
             else if(test == "pDiff"){
                  s<- sd(exp1)-sd(exp2);
                  theta <- 0;
                  for(k in 1:piter){
                     randExp<-sample(exp);
                     rand1Temp<-randExp[range1];
                     rand2Temp<-randExp[range2];
                     sdTemp1<-sd(rand1Temp);
                     sdTemp2<-sd(rand2Temp);
                     if(sdTemp1-sdTemp2>s){
                        theta <-theta +1;
                     }
                  }
                  statistic[g]<-s;
                  p.dv[g]<-theta /piter;
             }
              else if(test == "pMAD"){
                  s<- mad(exp1)-mad(exp2);
                  theta <- 0;
                  for(k in 1:piter){
                     randExp<-sample(exp);
                     rand1Temp<-randExp[range1];
                     rand2Temp<-randExp[range2];
                     madTemp1<-mad(rand1Temp);
                     madTemp2<-mad(rand2Temp);
                     if(madTemp1-madTemp2>s){
                        theta <-theta +1;
                     }
                  }
                  statistic[g]<-s;
                  p.dv[g]<-theta /piter;
             }
             else if(test == "pQn"){
                  s<- Qn(exp1)-Qn(exp2);
                  theta <- 0;
                  for(k in 1:piter){
                     randExp<-sample(exp);
                     rand1Temp<-randExp[range1];
                     rand2Temp<-randExp[range2];
                     QnTemp1<-Qn(rand1Temp);
                     QnTemp2<-Qn(rand2Temp);
                     if(QnTemp1-QnTemp2>s){
                        theta <-theta +1;
                     }
                  }
                  statistic[g]<-s;
                  p.dv[g]<-theta /piter;
             }
             else{
                print("ERROR: invalid option");
                return("ERROR: invalid option");
             }
             if(p.dv[g]!= -2 &&  p.dv[g]!= -3){
                  p.dv[g]<- 2*min((1-p.dv[g]), p.dv[g]);    #for two sided test
             } 
             #now, perform the t test
             tresult<-t.test(exp1, exp2, alternative="two.sided");
             p.de[g]<-tresult$p.value;
             t.value[g]<-tresult$statistic;
        }
   }
   result<-data.frame(actualDF1,actualDF2,statistic, t.value,p.dv,p.de);
   return(result);
}

###
# For Differential coexpression analysis
###

##
# calculate robust correlation coefficient using Pearon coefficient
#   x = an array of values, for sample 1
#   y = an array of values, for sample 2
#   outlier.tag = tha tag associated with outliers
##
robustCor<-function(x,y,outlier.tag=NA){
   tagged<-union(which(x==outlier.tag), which(y==outlier.tag));
   x<-x[setdiff(1:length(x), tagged)];
   y<-y[setdiff(1:length(y), tagged)];
   return(cor(x,y));
}

##
# calculate robust correlation coefficient using Pearon coefficient
#   data = a data matrix
#   method = the name of the coexpression measure, currently only "pearson" is implemented
#   outlier.rm = whether to remove outliers (TURE or FALSE)
#   r = the parameter for inter-quartile range filtering criteria
#   verbose = verbose?
##
robustCoexp<-function(data,method="pearson",outlier.rm=TRUE,r=3, verbose=TRUE){
   if(outlier.rm){
      numGene<-length(data[,1])
      for(i in 1:numGene){
         data[i,]<-detectOutlier(unlist(data[i,]),r);
      }
   }
   if(verbose){
      print("outlier detected");
   }
   size<-length(data[,1]);
   correlation<-array(-10,size*(size-1)/2);
   geneI<-array(-10,size*(size-1)/2);
   geneJ<-array(-10,size*(size-1)/2);
   k<-1;
   for(i in 1:(size-1)){
      for(j in (i+1):size){
         correlation[k]<-robustCor( t(data[i,]) , t(data[j,]), -1);
         geneI[k]<-i;
         geneJ[k]<-j;
         if(verbose){
            print(paste(geneI[k],geneJ[k]));
         }
         k<-k+1;
      }
   }
   result<-data.frame(geneI,geneJ,correlation);
   result

}

###
#  A method to simulate microarray data with realistic characteristics
#    g = number of genes
#    n = array of numbers indicating the number of arrays per group
#    type = type of changes
#    fc.dv = fold change in  variance ratio
#    fc.de = fold change in mean ratio
#    mu = mean of expression
#    sigma3 = variance of expression
#    distribution = type of distribution to be used
#    numOutlier = number of outliers
#    outlierConstant = the multiplier to the outlier values
###
simData<-function(g, n=c(), type=c("DV","DE"), fc.dv=1, fc.de=1, mu=7, sigma2=6, distribution="normal", numOutlier=0, outlierConstant=3){
   data<-c();
   mean<-mu;
   sd<-sqrt(sigma2);
   for(k in n){
      if(distribution=="normal"){
         data<-cbind(data, array(rnorm(g*k,mean,sd),c(g,k)));
      }
      else if(distribution=="uniform"){
         data<-cbind(data, array(runif(g*k,mean-sqrt(3)*sd,mean+sqrt(3)*sd),c(g,k)));
      }
      else if(distribution=="gamma"){
         data<-cbind(data, array(rgamma(g*k,shape=(mean^2/sd^2),rate=(mean/sd^2)),c(g,k)));
      }
      if(is.element("DV", type)){
         sd<-sd*fc.dv;
      }
      if(is.element("DE", type)){
         mean<-mean*fc.de;
      }
   }

   if(numOutlier>0){
      data[,runif(numOutlier, 1, sum(n))] <-  data[,1:numOutlier]*outlierConstant;
   }

   row.names(data)<- paste("gene", 1:g);
   return(data);
}


