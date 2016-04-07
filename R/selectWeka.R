selectWeka <- function(train,evaluator="CfsSubsetEval",search="BestFirst",n){
##

## evaluator:
## "CfsSubsetEval": Evaluates the worth of a subset of attributes by considering
##                  the individual predictive ability of each feature along with
##                  the degree of redundancy between them. Subsets of features
##                  that are highly correlated with the class while having low
##                  intercorrelation are preferred.
## "InfoGainAttributeEval":
## "SVMAttributeEval"

## search:
## "BestFirst": Searches the space of attribute subsets by greedy hillclimbing
##              augmented with a backtracking facility.
## "Ranker"

## n: the number of selected attributes.
##
  trainFile = tempfile("tempTrain")
  #require(foreign)
  write.arff(train,trainFile)
  tmpOut = tempfile("tempOut")
  if(search=="Ranker"){
      command = paste(.javaExecutable(), " -cp ", file.path(path.package("BioSeqClass"), "scripts", "weka.jar"),
                      " weka.attributeSelection.", evaluator, " -i ", trainFile,
                      " -s \"weka.attributeSelection.", search, " -N ", n, "\"", sep="" )
  }else{
      command = paste(.javaExecutable(), " -cp ", file.path(path.package("BioSeqClass"), "scripts", "weka.jar"),
                      " weka.attributeSelection.", evaluator, " -i ", trainFile,
                      " -s weka.attributeSelection.", search, sep="" )
  }
  tmp <- system(command,intern = TRUE)
  tmp <- tmp[grep("Selected attributes:",tmp)]
  tmp <- sub(" :.*","",sub("Selected attributes: ","",tmp))
  feature <- as.numeric(unlist(strsplit(tmp,split=",")))
  unlink(c(trainFile,tmpOut))
  feature
}
