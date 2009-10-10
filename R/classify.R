classify <- function(data,classifyMethod="libsvm",cv=10,
                           features, evaluator, search, n=200,
                           svm.kernel="linear",svm.scale=FALSE, 
                           svm.path, svm.options="-t 0",
                           knn.k=1,
                           nnet.size=2, nnet.rang=0.7, nnet.decay=0, nnet.maxit=100
                           ){
##
## data - a data frame. The last column is the class label and other columns are features.
##
## classifyMethod - a character for the method of classification.
##                  "libsvm": Support Vecotr Machine by LibSVM. Package "e1071" is required.
##                  "svmlight": Support Vecotr Machine by SVMLight. Package klaR is required.
##                  "NaiveBayes": Naive Bayes. Package klaR is required.
##                  "randomForest": random forest. Package randomForest is required.
##                  "knn": k Nearest Neighbor. Package class is required.
##                  "tree": Package tree is required.
##                  "nnet": neural net. Bundle VR is required.
##                  "rpart": Recursive Partitioning and Regression Trees. Package rpart is required.
##                  "ctree"
##                  "ctreelibsvm"
##                  "bagging"
##                  "svm-decisionTree"
##
## cv - Method for cross validation.
##      "5" : 
##      "10" :
##      "leave_one_out" : 
## 
## svm.kernel - a character for kernel function of SVM. 
##              "linear"
##              "polynomial"
##              "radial basis"
##              "sigmoid"
## svm.scale - A logical vector indicating the variables to be scaled.
## svm.path - a character for path to SVMlight binaries (required, if path is 
##            unknown by the OS).
## svm.options - Optional parameters to SVMlight.For further details see: 
##               ¡°How to use¡± on http://svmlight.joachims.org/. (e.g.: "-t 2 -g 0.1"))
## nnet.size - number of units in the hidden layer. Can be zero if there are 
##             skip-layer units. 
## nnet.rang - Initial random weights on [-rang, rang]. Value about 0.5 unless 
##             the inputs are large, in which case it should be chosen so that 
##             rang * max(|x|) is about 1. 
## nnet.decay - parameter for weight decay.
## nnet.maxit - maximum number of iterations.
##
  colnames(data)[ncol(data)] = "class"
  class = as.vector(unique(data[,ncol(data)]))
  if(cv=="leave_one_out"){
    cv = nrow(data)    
  }else{    
    count = floor(table(data[,ncol(data)])/cv)    
    index = lapply(class,function(y){which(data[,ncol(data)]==y)})
    names(index) = class    
  }
    
  result = list()
  for(i in 1:cv){
    if( cv==nrow(data) ){
      test = data[i,]
      train = data[-i,] 
    }else{
      tmp = unlist(lapply(class, function(x){index[[x]][(count[x]*(i-1)+1):(count[x]*i)]}) )
      test = data[tmp,]
      train = data[-tmp,]
    }
    
    if( missing(features) ){    
      if( missing(evaluator) ){
        features = 1:(ncol(train)-1)
      }else{
        features = selectWeka(train, evaluator, search, n)
      }
    }
    
    result[["features"]][[i]] = features
    train = train[,c(features,ncol(train))]
    test = test[,c(features,ncol(test))]
    
    if(classifyMethod=="knn"){
      testResult = classifyModelKNN(train, test[,-ncol(test)], knn.k)
    }else{
      if(classifyMethod=="ctreelibsvm"){
         testResult = classifyModelCTREELIBSVM(train, test[,-ncol(test)], svm.kernel, svm.scale)
      }else{
         model = switch(classifyMethod,
                        "libsvm" = classifyModelLIBSVM(train,svm.kernel,svm.scale),
                        "svmlight" = classifyModelSVMLIGHT(train,svm.path,svm.options="-t 0"),
                        "NaiveBayes" = classifyModelNB(train), 
                        "randomForest" = classifyModelRF(train),                 
                        "tree" = classifyModelTree(train),
                        "nnet" = classifyModelNNET(train,nnet.size,nnet.rang, 
                                 nnet.decay,nnet.maxit),
                        "rpart" = classifyModelRPART(train),
                        "ctree" = classifyModelCTREE(train),                       
                        "bagging" = classifyModelBAG(train),
                        stop(paste("classifyMethod",classifyMethod,"is not supported.
                        Has to be libsvm, svmlight, NaiveBayes, randomForest, 
                        tree, nnet, or rpart"))
                 )
         testResult = switch(classifyMethod,
                             "NaiveBayes" = predict(model, test[,-ncol(test)])$class,
                             "tree" = predict(model, test[,-ncol(test)],type="class"),
                             "nnet" = predict(model, test[,-ncol(test)],type="class"),
                             "rpart" = predict(model, test[,-ncol(test)],type="class"),
                             predict(model, test[,-ncol(test)])  
                      )
         }
    }
      
    if( cv==nrow(data) ){
      result[["test"]][[i]] = i
      result[["testPredict"]][[i]] = as.vector(testResult)
    }else{
      result[["test"]][[i]] = tmp
      perform = performance(testResult,test[,ncol(test)])  
      result[["performance"]][[i]] = perform
    }    
  }
  tmpName = names(result[["performance"]][[i]])
  tmpName1 = c("tp","tn","fp","fn")
  tmpName2 = setdiff(tmpName,tmpName1)
  if( cv==nrow(data) ){
    result[["totalPerformance"]] = performance(result[["testPredict"]],data[,ncol(data)])
  }else{
    result[["totalPerformance"]][tmpName1]=  
            apply(sapply(result[["performance"]],function(x){x[tmpName1]}), 
                  1, function(y){sum(y)}
            )
    result[["totalPerformance"]][tmpName2] = 
            apply(sapply(result[["performance"]],function(x){x[tmpName2]}), 
                  1, function(y){mean(y)}
            )
  }
  result[["featureNumber"]] = sort(table(unlist(
           result$features)),decreasing=T)         
  result
}
