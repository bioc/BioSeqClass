\name{performance}
\alias{performance}
\title{Performance Evaluation}
\description{
  Evaluate the performance of classification model. 
}
\usage{  
  performance(predictClass,factClass)
}
\arguments{
  \item{predictClass}{a factor of predicted classifications of training set, 
    comprising of "-1" or "+1".}
  \item{factClass}{a vector of true classifications of training set, 
    comprising of "-1" or "+1".}
}
\details{  
  \code{\link{performance}} evaluates the performance of classification model. It
    cacluates: tp (true positive), tn(ture negative), fp(false positive), 
    fn(false negative), prc(precision), sn(sensitivity), sp(specificity), 
    acc(accuracy), mcc(Matthews Correlation Coefficient), pc(Performance Coefficient).
}
\author{Hong Li}
\examples{
  ## read positive/negative sequence from files.
  tmpfile1 = file.path(path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  tmpfile2 = file.path(path.package("BioSeqClass"), "example", "acetylation_K.neg40.pep")
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
  ## Evaluate the performance of classification model
  performance(testClass,test[,ncol(test)])
}
