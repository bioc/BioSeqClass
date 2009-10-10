classifyModelLIBSVM <- function(train,svm.kernel="linear",svm.scale=FALSE ){
  #require(e1071)
  svm(train[,-ncol(train)],factor(train[,ncol(train)]),
      kernel=svm.kernel,scale=svm.scale)
}

classifyModelSVMLIGHT <- function(train,svm.path,svm.options="-t 0" ){
  if( missing(svm.path) ){
    stop("Parameter 'svm.path' must not be missing")
  }
  #require(klaR)
  svmlight(train[,-ncol(train)],factor(train[,ncol(train)]),
           pathsvm=svm.path, svm.options=svm.options)
}

classifyModelNB <- function(train){
  #require(klaR)  
  model <- NaiveBayes(class ~., data.frame(train))
}

classifyModelRF <- function(train){
  #require(randomForest)  
  model <- randomForest(train[,-ncol(train)],factor(train[,ncol(train)]) )
}

classifyModelKNN <- function(train, test, knn.k=1){
  #require(class)
  if( is.null(nrow(test)) ){
    predictResult <- knn(train[,-ncol(train)], matrix(test), 
      factor(train[,ncol(train)]), k=knn.k)
  }else{
    predictResult <- knn(train[,-ncol(train)], test, 
      factor(train[,ncol(train)]), k=knn.k)
  }
  predictResult
}

classifyModelTree <- function(train){
  #require(tree)
  model <- tree(class ~., data.frame(train))
}

classifyModelNNET <- function(train, nnet.size=2, nnet.rang=0.7, nnet.decay=0, 
                     nnet.maxit=100){
  #require(nnet)
  model <- nnet(class ~., train, size=nnet.size, rang=nnet.rang, 
           decay=nnet.decay, maxit=nnet.maxit)
}

classifyModelRPART <- function(train){
  #require(rpart)
  model <- rpart(class ~., train, method="class")
}
  
classifyModelCTREE <- function(train){ 
  #require(party)
  model <- ctree(class ~., train)
}

classifyModelCTREELIBSVM<- function(train, test, svm.kernel="linear",svm.scale=FALSE){ 
  #require(party)
  #require(e1071)
  treeModel <- ctree(class ~., train)
  pdf(paste(tempfile("CTREE"),"pdf",sep="."))
  plot(treeModel)
  dev.off()
  node <- where(treeModel)  
  multiModel <- tapply(1:length(node),node, function(x){                
                  svm(train[x,-ncol(train)],factor(train[x,ncol(train)]),
                  kernel=svm.kernel,scale=svm.scale)
                })                   
  testGroup <- tapply(1:nrow(test), where(treeModel,test), function(x){
                 test[x,]
               }) 
  predictResult <- unlist(sapply(names(testGroup), function(x){
                     tmp <- testGroup[[x]]
                     predict(multiModel[x], tmp[,-ncol(tmp)])
                   }))
  names(predictResult)<- unlist(sapply(names(testGroup), function(x){ 
                         rownames(testGroup[[x]]) }))
  predictResult[rownames(test)]                       
}

classifyModelBAG<- function(train){
  #require(ipred)  
  model <- bagging(class ~., train)
}

