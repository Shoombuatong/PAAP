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

customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes
######### Optimized parameter
control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(1:10), .ntree=c(100,200,300,400,500))
custom <- train(Class~., data=training , method=customRF, metric=c("Accuracy"), tuneGrid=tunegrid, trControl=control)

Model = randomForest(Class ~ ., training, ntree= 300,mtry =)

#######Feature extraction
xtest <- read.fasta('AHTP203.fasta', seqtype="AA", as.string = TRUE)###read data
aactest<- t(sapply(xtest, extractAAC))
dpctest <- t(sapply(xtest, extractDC))
data <- data.frame(aactest,dpctest)

#######Predicting unknown sequences
data.frame(Prediction= predict(Model,data))
