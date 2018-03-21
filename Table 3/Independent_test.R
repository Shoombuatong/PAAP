#######set directory
setwd('D:\\Peptide prediction\\Antihypertensive peptides\\CLassification\\Backup')
#######Load package
library(caret)
library(randomForest)
library(protr)
library(seqinr)
library(Interpol)

#######Building RF model
A <- read.fasta('AHTP2132.fasta', seqtype="AA", as.string = TRUE)
D = read.csv("Class2132.csv", header = TRUE) 
aac<- t(sapply(A, extractAAC))
dpc <- t(sapply(A, extractDC))
training = data.frame(aac,dpc,Class = D[,ncol(D)])
Model = randomForest(Class ~ ., training, ntree= 300,mtry =1)

#######Feature extraction
xtest <- read.fasta('AHTP203.fasta', seqtype="AA", as.string = TRUE)###read data
aactest<- t(sapply(xtest, extractAAC))
dpctest <- t(sapply(xtest, extractDC))
data <- data.frame(aactest,dpctest)

#######Predicting unknown sequences
data.frame(Prediction= predict(Model,data))
