\name{classify}
\alias{classify}
\title{Classification}
\description{
  Classification with cross validation and selected features
}
\usage{
  classify(data,classifyMethod="libsvm",cv=10,
                 features, evaluator, search, n=200,        
                 svm.kernel="linear",svm.scale=FALSE, 
                 svm.path, svm.options="-t 0",
                 knn.k=1,
                 nnet.size=2, nnet.rang=0.7, nnet.decay=0, nnet.maxit=100)
}
\arguments{
  \item{data}{a data frame including the feature matrix and class label. The last
    column is a vector of class label comprising of "-1" or "+1"; 
    Other columns are features.} 
  \item{classifyMethod}{a string for the classification method. This must be one 
    of the strings "libsvm", "svmlight", "NaiveBayes", "randomForest", "knn",
    "tree", "nnet", "rpart", "ctree", "ctreelibsvm", "bagging".} 
  \item{cv}{an integer for the time of cross validation, or a string "leave\_one\_out" 
    for the jacknife test.} 
  \item{features}{an integer vector for the index of interested columns in data,
    which will be used as features for build classification model.}
  \item{evaluator}{a string for the feature selection method used by WEKA. This 
    must be one of the strings "CfsSubsetEval", "ChiSquaredAttributeEval", 
    "InfoGainAttributeEval", or "SVMAttributeEval".} 
  \item{search}{a string for the search method used by WEKA. This must be one 
    of the strings "BestFirst" or "Ranker".} 
  \item{n}{an integer for the number of selected features.} 
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
}
\details{  
  \code{\link{classify}} employ feature selction method in Weka and 
  diverse classification model in other R packages to perfrom classification.
  "Cross Validation" is controlled by parameter "cv"; 
  "Feature Selection" is controlled by parameter "features", "evaluator", "search", 
  and "n";
  "Classification Model Building" is controlled by parameter "classifyMethod". 

  Parameter "evaluator" supportes three feature selection methods provided by WEKA:
  "CfsSubsetEval": Evaluate the worth of a subset of attributes by considering 
                   the individual predictive ability of each feature along with 
                   the degree of redundancy between them.
  "ChiSquaredAttributeEval": Evaluate the worth of an attribute by computing the 
                             value of the chi-squared statistic with respect to 
                             the class.
  "InfoGainAttributeEval": Evaluate attributes individually by measuring 
                           information gain with respect to the class.
  "SVMAttributeEval": Evaluate the worth of an attribute by using an SVM classifier. 
                      Attributes are ranked by the square of the weight assigned 
                      by the SVM. Attribute selection for multiclass problems is 
                      handled by ranking attributes for each class seperately 
                      using a one-vs-all method and then "dealing" from the top 
                      of each pile to give a final ranking.
                      
  Parameter "search" supportes three feature subset search methods provided by WEKA:
  "BestFirst":  Searches the space of attribute subsets by greedy hillclimbing 
                augmented with a backtracking facility. Setting the number of 
                consecutive non-improving nodes allowed controls the level of 
                backtracking done. Best first may start with the empty set of 
                attributes and search forward, or start with the full set of 
                attributes and search backward, or start at any point and search 
                in both directions (by considering all possible single attribute 
                additions and deletions at a given point).
  "Ranker": Ranks attributes by their individual evaluations.

  Parameter "classifyMethod" supports multiple classification model:
  "libsvm": Employ \code{\link{classifyModelLIBSVM}} to perform Support Vecotr 
            Machine by LibSVM. Package "e1071" is required.
  "svmlight": Employ \code{\link{classifyModelSVMLIGHT}} to Support Vecotr Machine 
              by SVMLight. Package "klaR" is required.
  "NaiveBayes": Employ \code{\link{classifyModelNB}} to perform Naive Bayes 
                classification. Package "klaR" is required.
  "randomForest": Employ \code{\link{classifyModelRF}} to perform random forest
                  classification. Package "randomForest" is required.
  "knn": Employ \code{\link{classifyModelKNN}} to perform k Nearest Neighbor 
         algorithm. Package "class" is required.
  "tree": Employ \code{\link{classifyModelTree}} to perform tree classification.
          Package "tree" is required.
  "nnet": Employ \code{\link{classifyModelNNET}} to perform neural net algorithm. 
          Bundle "VR" is required.
  "rpart": Employ \code{\link{classifyModelRPART}} to perform Recursive 
          Partitioning and Regression Trees. Package "rpart" is required.
  "ctree": Employ \code{\link{classifyModelCTREE}} to perform Conditional 
           Inference Trees. Package "party" is required.
  "ctreelibsvm": Employ \code{\link{classifyModelCTREELIBSVM}} to combine Conditional 
                 Inference Trees and Support Vecotr Machine for classification. 
                 For each node in the tree, one SVM model will be constructed using
                 train data in this node. Test data will be firstly classified 
                 to one node of the tree, and then use corresponding SVM to do 
                 classification. Package "party" and "e1071" is required.
  "bagging": Employ \code{\link{classifyModelBAG}} to perform bagging for 
             classification trees. Package "ipred" is required.
             
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
  
  ## Use LibSVM and 5-cross-validation to classify.
  LIBSVM_CV5 = classify(data,classifyMethod="libsvm",cv=5,
                 svm.kernel="linear",svm.scale=FALSE)
  ## Features selection is done by envoking "CfsSubsetEval" method in WEKA.               
  FS_LIBSVM_CV5 = classify(data,classifyMethod="libsvm",cv=5,evaluator="CfsSubsetEval",
                 search="BestFirst",svm.kernel="linear",svm.scale=FALSE)    
  
if(interactive()){

  KNN_CV5 = classify(data,classifyMethod="knn",cv=5,knn.k=1)  
  
  RF_CV5 = classify(data,classifyMethod="randomForest",cv=5)
  
  TREE_CV5 = classify(data,classifyMethod="tree",cv=5)
  
  NNET_CV5 = classify(data,classifyMethod="nnet",cv=5)
  
  RPART_CV5 = classify(data,classifyMethod="rpart",cv=5,evaluator="")
  
  CTREE_CV5 = classify(data,classifyMethod="ctree",cv=5,evaluator="")  
  
  BAG_CV5 = classify(data,classifyMethod="bagging",cv=5,evaluator="")  
       
}
}
  