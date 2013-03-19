\name{featureEvaluate}
\alias{featureEvaluate}
\title{Evaluate Different Feature Coding Schemas}
\description{
  Feature sets from different feature coding schemas are used as input of classification 
  models, and the model performance are given in the result.
}
\usage{
  featureEvaluate(seq, classLable, fileName, ele.type, featureMethod, 
            cv=10, classifyMethod="libsvm",             
            group=c("aaH", "aaV", "aaZ", "aaP", "aaF", "aaS", "aaE"), k, g, 
            hydro.methods=c("kpm", "SARAH1"), hydro.indexs=c("hydroE", "hydroF", "hydroC"), 
            aaindex.name, n, d, w=0.05, start.pos, stop.pos, psiblast.path, 
            database.path, hmmpfam.path, pfam.path, Evalue=10^-5,       
            na.type="all", na.strand="all", diprodb.method="all", diprodb.type="all",     
            svm.kernel="linear", svm.scale=FALSE, svm.path, svm.options="-t 0",
            knn.k=1, nnet.size=2, nnet.rang=0.7, nnet.decay=0, nnet.maxit=100)
}
\arguments{
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.}
  \item{classLable}{a factor or vector for the class lable of sequences in seq.}
  \item{fileName}{a string for the output file name.}  
  \item{ele.type}{a string for the type of biological sequence. This must be 
    one of the strings "rnaBase", "dnaBase", "aminoacid" or "aminoacid2".} 
  \item{featureMethod}{a string vector for the name of feature coding. The 
    alternative names are "Binary", "CTD", "FragmentComposition", 
    "GapPairComposition", "CKSAAP", "Hydro", "ACH", "AAindex", "ACI", 
    "ACF", "PseudoAAComp", "PSSM", "DOMAIN", "BDNAVIDEO", and "DIPRODB".} 
  \item{classifyMethod}{a string for the classification method. This must be one 
    of the strings "libsvm", "svmlight", "NaiveBayes", "randomForest", "knn",
    "tree", "nnet", "rpart", "ctree", "ctreelibsvm", "bagging".} 
  \item{cv}{an integer for the time of cross validation, or a string 
    "leave\_one\_out" for the jacknife test.} 
  \item{group}{a string vector for the group of amino acids. This alternative 
    groups are: "aaH", "aaV", "aaZ", "aaP", "aaF", "aaS" or "aaE".}
  \item{k}{an integer indicating the length of sequence fragment (k>=1).} 
  \item{g}{an integer indicating the distance between two aminoacids/bases
    (g>=0).}  
  \item{hydro.methods}{a string vector for the methods of coding protein 
    hydrophobic effect. This alternative groups are: "kpm" or "SARAH1".}
  \item{hydro.indexs}{a string vector for the methods of coding protein 
    hydrophobic effect. This alternative groups are: "hydroE", "hydroF" or "hydroC".}
  \item{aaindex.name}{a string for the name of physicochemical and biochemical 
    properties in AAindx.}    
  \item{n}{an integer used as paramter of \code{\link{featureACF}} (1<=n<=L-2, 
    L is the the length of sequence). featureACF takes the auto-correlation 
    between fragment X(1)...X(L-m) and X(m+1)...X(L) (1<=m<=n) as features. }
  \item{d}{an integer used as paramter of \code{\link{featurePseudoAAComp}} 
    (d>=1). Coupling between amino acids X(i) and X(i+d) are considered as features.}
  \item{w}{a numeric value for the weight factor of sequence order effect in 
    \code{\link{featurePseudoAAComp}}. }
  \item{start.pos}{a integer vector denoting the start position of the fragment window.
    If it is missing, it is 1 by default.}
  \item{stop.pos}{a integer vector denoting the stop position of the fragment window.
    If it is missing, it is the length of sequence by default.} 
  \item{psiblast.path}{a string for the path of PSI-BLAST program blastpgp. blastpgp
    will be employed to iteratively search database and generate position-specific 
    scores for each position in the alignment.}
  \item{database.path}{a string for the path of formatted protein database. Database
    can be formatted by formatdb program.}
  \item{hmmpfam.path}{a string for the path of hammpfam program in HMMER. 
    hammpfam will be employed to predict domains using models in Pfam database.}
  \item{pfam.path}{a string for the path of pfam domain database.}
  \item{Evalue}{a numeric value for the E-value cutoff of perdicted Pfam domain.}
  \item{na.type}{a string for nucleic acid type. It must be "DNA", "DNA/RNA", "RNA", 
    or "all".}
  \item{na.strand}{a string for strand information. It must be "double", "single", 
    or "all".}
  \item{diprodb.method}{a string for mode of property determination. It can be 
    "experimental", "calculated", or "all".}
  \item{diprodb.type}{a string for property type. It can be "physicochemical", 
    "conformational", "letter based", or "all".}
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
  \code{\link{featureEvaluate}} can test feature coding methods for short 
  peptide, protein, DNA or RNA. 
  It returns a ranked list based on the accuracy of classification result. 
  Each element in the list has three components: "data", "model", and "performance".
  "data" is a data.frame object, which stores feature matrix and its last column 
  is the class label. "model" is a vector for feature coding method, which 
  contains 6 elements: "Feature\_Function", "Feature\_Parameter", 
  "Feature\_Number", "Model", "Model\_Parameter", and "Cross_Validataion". 
  "performance" is a vector for the performance result of classification model, 
  which contains 10 elements: "tp", "tn", "fp", "fn", "prcc", "sn", "sp", "acc", 
  "mcc", "pc".   
}
\author{Hong Li}
\examples{
  ## read positive/negative sequence from files.
  tmpfile1 = file.path(path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  tmpfile2 = file.path(path.package("BioSeqClass"), "example", "acetylation_K.neg40.pep")
  posSeq = as.matrix(read.csv(tmpfile1,header=FALSE,sep="\t",row.names=1))[,1]
  negSeq = as.matrix(read.csv(tmpfile2,header=FALSE,sep="\t",row.names=1))[,1]
  seq=c(posSeq,negSeq)
  classLable=c(rep("+1",length(posSeq)),rep("-1",length(negSeq)) ) 
  if(interactive()){    
   ## test various feature coding methods.
   ## it may be time consuming.
    fileName = tempfile()
    testFeatureSet = featureEvaluate(seq, classLable, fileName, ele.type="aminoacid", 
              featureMethod=c("Binary", "CTD", "FragmentComposition", "GapPairComposition", 
              "Hydro"), cv=5, classifyMethod="libsvm",             
              group=c("aaH", "aaV", "aaZ", "aaP", "aaF", "aaS", "aaE"), k=3, g=7, 
              hydro.methods=c("kpm", "SARAH1"), hydro.indexs=c("hydroE", "hydroF", "hydroC") )            
    summary = read.csv(fileName,sep="\t",header=T)
    fix(summary)    
    
    ## Evaluate features from different feature coding functions
    feature.index = 1:5
    tmp <- testFeatureSet[[1]]$data
    colnames(tmp) <- paste(testFeatureSet[[feature.index[1]]]$model["Feature_Function"],testFeatureSet[[feature.index[1]]]$model["Feature_Parameter"],colnames(tmp),sep=" ; ")      
    data <- tmp[,-ncol(tmp)]  
    for(i in 2:length(feature.index) ){
      tmp <- testFeatureSet[[feature.index[i]]]$data
      colnames(tmp) <- paste(testFeatureSet[[feature.index[i]]]$model["Feature_Function"],testFeatureSet[[feature.index[i]]]$model["Feature_Parameter"],colnames(tmp),sep=" ; ")    
      data <- data.frame(data, tmp[,-ncol(tmp)] )
    }
    name <- colnames(data)
    data <- data.frame(data, tmp[,ncol(tmp)] )
    ## feature forward selection by 'cv_FFS_classify'
    ## it is very time consuming.
    combineFeatureResult = fsFFS(data,stop.n=50,classifyMethod="knn",cv=5)  
    tmp = sapply(combineFeatureResult,function(x){c(length(x$features),x$performance["acc"])})
    plot(tmp[1,],tmp[2,],xlab="featureNumber",ylab="Accuracy",main="result of FFS_KNN",pch=19)
    lines(tmp[1,],tmp[2,])
    
    ## compare the prediction accuracy based on different feature coding methods and different classification models.
    ## it is very time consuming.
    testResult = lapply(c("libsvm", "randomForest", "knn", "tree"), 
      function(x){
              tmp = featureEvaluate(seq, classLable, fileName = tempfile(), 
              ele.type="aminoacid", featureMethod=c("Binary", "CTD", "FragmentComposition", 
              "GapPairComposition", "Hydro"), cv=5, classifyMethod=x,             
              group=c("aaH", "aaV", "aaZ", "aaP", "aaF", "aaS", "aaE"), k=3, g=7, 
              hydro.methods=c("kpm", "SARAH1"), hydro.indexs=c("hydroE", "hydroF", "hydroC") );
              sapply(tmp,function(y){c(y$model[["Feature_Function"]], y$model[["Feature_Parameter"]], y$model[["Model"]], y$performance[["acc"]])})
    })
    tmpFeature = as.factor(c(sapply(testResult,function(x){apply(x[1:2,],2,function(y){paste(y,collapse="; ")})})))
    tmpModel = as.factor(c(sapply(testResult,function(x){x[3,]})))
    tmp1 = data.frame(as.integer(tmpFeature), as.integer(tmpModel), as.numeric(c(sapply(testResult,function(x){x[4,]}))) )
    require(scatterplot3d)
    s3d=scatterplot3d(tmp1,color=c("red","blue","green","yellow")[tmp1[,2]],pch=19,
        xlab="Feature Coding", ylab="Classification Model", 
        zlab="Accuracy under 5-fold cross validation",lab=c(10,6,7),     
        y.ticklabs=c("",as.character(sort(unique(tmpModel))),"") )
  }
}
  
