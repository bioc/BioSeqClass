BioSeqClass <- function(version="1.5.0"){    
  
  depends = c(
         "ipred","e1071","klaR","randomForest","class","tree","nnet",
         "rpart","party","foreign","Biobase","Biostrings")  
  tmp <- sapply(depends,function(x){ 
    if (!require(x,quietly=TRUE,character.only=TRUE) ){
      print(paste("Installing depended package",x,"..."))
      if( any(x==c("ipred","e1071","klaR","randomForest","class","tree","nnet",
         "rpart","party","foreign")) ){
        install.packages(x,repos="http://stat.ethz.ch/CRAN",quiet=T)
      }
      if( any(x==c("Biobase","Biostrings")) ){
        if (!requireNamespace("BiocManager", quietly=TRUE))
            install.packages("BiocManager")
        BiocManager::install(x)
      }
      require(x,quietly=TRUE,character.only=TRUE) || stop(paste("Package",x,
      "is not correctly installed"))   
    }          
  }) 
      
  OS <- .Platform$OS.type
  if(OS == "windows"){
    tmpFile <- paste("http://www.biosino.org/download/BioSeqClass/BioSeqClass_",version,".zip", sep="")
    tmpPath <- file.path(tempdir(),paste("BioSeqClass_",version,".zip",sep=""))
  }else{
    tmpFile <- paste("http://www.biosino.org/download/BioSeqClass/BioSeqClass_",version,".tar.gz", sep="")
    tmpPath <- file.path(tempdir(),paste("BioSeqClass_",version,".tar.gz",sep=""))
  }    
  download.file(tmpFile, tmpPath)
  install.packages(tmpPath,repos = NULL) 
}

