\name{selectFFS}
\alias{selectFFS}
\title{Feature Selection by FFS}
\description{
  feature forward selection.
}
\usage{
  selectFFS(data, accCutoff, stop.n,
            classifyMethod="knn",cv=10)  
}
\arguments{
  \item{data}{a data frame including the feature matrix and class label. The last
    column is a vector of class label comprising of "-1" or "+1"; 
    Other columns are features.}
  \item{accCutoff}{a numeric indicating the minimum difference of accuracy 
    between two models in \code{\link{selectFFS}}. Feature subsets will stop
    increasing when the difference of accuracy is samll than accCutoff.} 
  \item{stop.n}{number of selected features by \code{\link{selectFFS}}.}   
  \item{classifyMethod}{a string for the classification method. This must be one 
    of the strings "libsvm", "svmlight", "NaiveBayes", "randomForest", "knn",
    "tree", "nnet", "rpart", "ctree", "ctreelibsvm", "bagging".} 
  \item{cv}{an integer for the time of cross validation, or a string "leave\_one\_out" 
    for the jacknife test.}   
}
\details{  
  \code{\link{selectFFS}} uses FFS (Feature Forword Selection) method to 
  increase feature, and use NNA (Neareast Neighbor Analysis) to evaluate 
  the performance of feature subset. Two conditions are used to stop feature 
  increasing: control the difference of accuracy between two models; control 
  the number of selected features by Parameter "stop.n".
}
\author{Hong Li}
\examples{
  ## read positive/negative sequence from files.
  tmpfile1 = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  tmpfile2 = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.neg40.pep")
  posSeq = as.matrix(read.csv(tmpfile1,header=FALSE,sep="\t",row.names=1))[,1]
  negSeq = as.matrix(read.csv(tmpfile2,header=FALSE,sep="\t",row.names=1))[,1]
  seq=c(posSeq,negSeq)
  classLable=c(rep("+1",length(posSeq)),rep("-1",length(negSeq)) ) 
  data = data.frame(featureBinary(seq),classLable)
  
if(interactive()){  
  ## Use KNN to evaluate the performance of feature subset, 
  ## and use Feature Forword Selection method to increase feature.
  # If the difference of accuracy between two models is less than 0.01, feature 
  # selection will stop.
  FFS_NNA_CV5 = selectFFS(data,accCutoff=0.01,classifyMethod="knn",cv=5)
  # If 20 features have been selected, feature selection will stop.
  FFS_NNA_CV5 = selectFFS(data,stop.n=3,classifyMethod="knn",cv=5)
  # If any one condiction is satisfied, feature selection will stop.
  FFS_NNA_CV5 = selectFFS(data,accCutoff=0.001,stop.n=100,classifyMethod="knn",cv=5)   
}
}
  