#######set directory
setwd('D:\\Peptide prediction\\Antihypertensive peptides\\CLassification')
#######Load package
library(caret)
library(randomForest)
library(protr)
library(seqinr)
library(Interpol)

##### Read data
x <- read.fasta('pentaall.fasta', seqtype="AA", as.string = TRUE)
D = read.csv("class_pentapeptide.csv", header = TRUE) 

###### Feature extraction
m = length(x)
aac <- t(sapply(x, extractAAC))
internal = data.frame(aac,Class = D[,ncol(D)])

ind= c(2,3,5,7,9,11,13,15,17,20)
n = ncol(internal)-1
gini = matrix(nrow = n, ncol = 10)
meangini = matrix(nrow = 400, ncol = 1)

for (i in 1:10){
RF<-randomForest(Class~.,data=internal,ntree=100,mtry=ind[i],importance=TRUE)
gini[,i] = RF$ importance[,4]
}

for (i in 1:400){
meangini[i,] = mean(gini[i,])
}
