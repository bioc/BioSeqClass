\name{model}
\alias{model}
\alias{classifyModelLIBSVM}
\alias{classifyModelSVMLIGHT}
\alias{classifyModelNB}
\alias{classifyModelRF}
\alias{classifyModelKNN}
\alias{classifyModelTree}
\alias{classifyModelRPART}
\alias{classifyModelCTREE}
\alias{classifyModelCTREELIBSVM}
\alias{classifyModelBAG}
\alias{classifyModelNNET}
\title{Classification Models}
\description{
  These functions build various classification models. 
}
\usage{  
  classifyModelLIBSVM(train,svm.kernel="linear",svm.scale=FALSE)
  classifyModelSVMLIGHT(train,svm.path,svm.options="-t 0")
  classifyModelNB(train)
  classifyModelRF(train)
  classifyModelKNN(train, test, knn.k=1)
  classifyModelTree(train)
  classifyModelNNET(train, nnet.size=2, nnet.rang=0.7, nnet.decay=0, nnet.maxit=100)
  classifyModelRPART(train)
  classifyModelCTREE(train)
  classifyModelCTREELIBSVM(train, test, svm.kernel="linear",svm.scale=FALSE)
  classifyModelBAG(train)  
}
\arguments{
  \item{train}{a data frame including the feature matrix and class label.
    The last column is a vector of class label comprising of "-1" or "+1"; 
    Other columns are features.}
  \item{svm.kernel}{a string for kernel function of SVM.} 
  \item{svm.scale}{a logical vector indicating the variables to be scaled.} 
  \item{svm.path}{a character for path to SVMlight binaries (required, if path 
    is unknown by the OS).} 
  \item{svm.options}{Optional parameters to SVMlight. For further details see: 
    "How to use" on http://svmlight.joachims.org/. (e.g.: "-t 2 -g 0.1")) } 
  \item{nnet.size}{number of units in the hidden layer. Can be zero if there are 
    skip-layer units.} 
  \item{nnet.rang}{Initial random weights on [-rang, rang]. Value about 0.5 unless 
    the inputs are large, in which case it should be chosen so that 
    rang * max(|x|) is about 1. } 
  \item{nnet.decay}{parameter for weight decay.} 
  \item{nnet.maxit}{maximum number of iterations.}   
  \item{knn.k}{number of neighbours considered in function \code{\link{classifyModelKNN}}.}  
  \item{test}{a data frame including the feature matrix and class label.
    The last column is a vector of class label comprising of "-1" or "+1"; 
    Other columns are features.}  
}
\details{  
  \code{\link{classifyModelLIBSVM}} builds support vector machine model by LibSVM.
    R package "e1071" is needed.
 
  \code{\link{classifyModelSVMLIGHT}} builds support vector machine model by SVMlight.  
    R package "klaR" is needed.
    
  \code{\link{classifyModelNB}} builds naive bayes model. R package "klaR" is 
    needed.

  \code{\link{classifyModelRF}} builds random forest model. R package 
    "randomForest" is needed.
    
  \code{\link{classifyModelKNN}} builds k-nearest neighbor model. R package 
    "class" is needed.
    
  \code{\link{classifyModelTree}} builds tree model. R package "class" is needed.

  \code{\link{classifyModelRPART}} builds recursive partitioning trees model. R 
    package "rpart" is needed.
   
  \code{\link{classifyModelCTREE}} builds conditional inference trees model. R 
    package "party" is needed.
   
  \code{\link{classifyModelCTREELIBSVM}} combines conditional inference trees and 
    support vecotr machine. R package "party" and "e1071" is needed.
    
  \code{\link{classifyModelBAG}} uses bagging method. R package "ipred" is needed.  
}
\author{Hong Li}
\examples{
  ## read positive/negative sequence from files.
  tmpfile1 = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  tmpfile2 = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.neg40.pep")
  posSeq = as.matrix(read.csv(tmpfile1,header=FALSE,sep="\t",row.names=1))[,1]
  negSeq = as.matrix(read.csv(tmpfile2,header=FALSE,sep="\t",row.names=1))[,1]
  data = data.frame(rbind(featureBinary(posSeq,elements("aminoacid")), 
       featureBinary(negSeq,elements("aminoacid")) ),
       class=c(rep("+1",length(posSeq)),
       rep("-1",length(negSeq))) )
  
  ## sample train and test data
  tmp = c(sample(1:length(posSeq),length(posSeq)*0.8), 
    sample(length(posSeq)+(1:length(negSeq)),length(negSeq)*0.8))
  train = data[tmp,]
  test = data[-tmp,]
  
  ## Build classification model using training data
  model1 = classifyModelLIBSVM(train,svm.kernel="linear",svm.scale=FALSE)
  ## Predict test data by classification model
  testClass = predict(model1, test[,-ncol(test)])  
}
