#######set directory
setwd('D:\\Peptide prediction\\Antihypertensive peptides\\CLassification\\Backup')
#######Load package
library(caret)
library(randomForest)
library(protr)
library(seqinr)
library(Interpol)
library(pls)

##### Set parameter for PseAAC
para1 = 2
para2 = 0.1
##### Read data
x <- read.fasta('pentaall.fasta', seqtype="AA", as.string = TRUE)
D = read.csv("class_pentapeptide.csv", header = TRUE) 

###### Feature extraction
m = length(x)
aac <- t(sapply(x, extractAAC))
dpc <- t(sapply(x, extractDC))
paac <- matrix(nrow = m, ncol = 20 + para1)
for(i in 1:m){ 
paac[i, ]= extractPAAC(x[[i]][1],lambda = para1 , w = para2 , props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

######## Select feature for constructing predictor 
internal = data.frame(aac,dpc,Class = D[,ncol(D)])

######## The customRF is used for optimizing lamda abd weight
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
control <- trainControl(method="cv", number=5)
tunegrid <- expand.grid(.mtry=c(1:10), .ntree=c(100,200,300,400,500))
custom <- train(Class~., data=internal , method=customRF, metric=c("Accuracy"), tuneGrid=tunegrid, trControl=control)

######Loop for LOOCV
k <- nrow(internal);
Resultcv <- 0;
folds <- cvsegments(nrow(internal), k);

for (fold in 1:k){
  currentFold <- folds[fold][[1]];
  RF = randomForest(Class ~ ., internal[-currentFold,], ntree= as.numeric(custom$ bestTune[2]),mtry = as.numeric(custom$ bestTune[1]),
  orm.votes=TRUE,keep.forest=TRUE, importance=TRUE) ## Building RF model
  pred = predict(RF, internal[currentFold,])
  Resultcv <- Resultcv + table(true=internal[currentFold,]$Class, pred=pred);   

################### Performance report
data = Resultcv
	ACCtr = (data[1]+data[4])/(data[1]+data[2]+data[3]+data[4])*100
	SENStr  =  (data[1]/(data[1]+data[2]))*100
	SPECtr = (data[4])/(data[3]+data[4])*100
	MCC1      = (data[1]*data[4]) - (data[2]*data[3])
	MCC2      =  (data[4]+data[2])*(data[4]+data[3])
	MCC3      =  (data[1]+data[2])*(data[1]+data[3])
	MCC4	=  sqrt(MCC2)*sqrt(MCC3)
	MCCtr  = MCC1/MCC4
}
result = data.frame (ACCtr,SENStr,SPECtr,MCCtr)

