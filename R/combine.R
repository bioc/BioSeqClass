combine <- function(seq, classLable, fileName, ele.type, featureMethod, 
            cv=10, classifyMethod="libsvm",             
            group=c("aaH", "aaV", "aaZ", "aaP", "aaF", "aaS", "aaE"), k, g, 
            hydro.methods=c("kpm", "SARAH1"), hydro.indexs=c("hydroE", "hydroF", "hydroC"), 
            aaindex.name, n, d, w=0.05, start.pos, stop.pos, psiblast.path, 
            database.path, hmmpfam.path, pfam.path, Evalue=10^-5,       
            na.type="all", na.strand="all", diprodb.method="all", diprodb.type="all",     
            svm.kernel="linear", svm.scale=FALSE, svm.path, svm.options="-t 0",
            knn.k=1, nnet.size=2, nnet.rang=0.7, nnet.decay=0, nnet.maxit=100
){  
  if( any(c(missing(seq), missing(classLable), missing(fileName), 
    missing(ele.type), missing(featureMethod),
    is.null(seq), is.null(classLable), is.null(fileName), 
    is.null(ele.type), is.null(featureMethod) ))  ){
     stop("Parameter 'seq', 'classLable', 'fileName', 'ele.type' and 
     'featureMethod' must not be missing")
  }
  if( !any(ele.type==c("rnaBase","dnaBase","aminoacid","aminoacid2")) ){
     stop("parameter 'ele.type' has to be: rnaBase, dnaBase, aminoacid 
             or aminoacid2")
  }
  if( length(setdiff(featureMethod, c("Binary", "CTD", "FragmentComposition", 
    "GapPairComposition", "CKSAAP", "Hydro", "ACH", "AAindex", "ACI", 
    "ACF", "PseudoAAComp", "PSSM", "DOMAIN", "BDNAVIDEO", "DIPRODB")))>0 ){
     stop("parameter 'featureMethod' has to be: ")
  }
  if( ele.type!="aminoacid" & length(intersect(featureMethod, c("Hydro", "ACH", "AAindex", 
    "ACI", "ACF", "PseudoAAComp", "PSSM", "DOMAIN")))>0 ){
     stop("Only 'ele.type=aminoacid' can use feature coding methods: 
     Hydro, ACH, AAindex, ACI, ACF, PseudoAAComp, PSSM and DOMAIN")     
  }  
  L = unique(sapply(seq,function(x){length(unlist(strsplit(x,split="")))}))  
  if( length(L)>1 ){
    if( length(intersect(featureMethod,c("Binary","Hydro","ACH",
      "AAindex")))>0 ){
      stop("Feature Coding Method (", paste(c("Binary","Hydro","ACH",
      "AAindex"), collapse="; ") ,") needs input sequences have equal length")
    }
  }
  if( sum(L%%2==0)>0 & any(featureMethod=="ACI") ){
    stop("Feature Coding Method 'ACI' needs sequences with odd elements")
  }
  if( !any(classifyMethod==c("libsvm", "svmlight", "NaiveBayes", "randomForest",
     "knn", "tree", "nnet", "rpart", "ctree", "ctreelibsvm", "bagging")) ){
     stop("parameter 'classifyMethod' has to be: libsvm, svmlight,      
       NaiveBayes, randomForest, knn, tree, nnet, rpart, ctree, ctreelibsvm, 
       or bagging")          
  }
  if( classifyMethod=="svmlight" ){
    if( missing(svm.path) ){
      stop("Parameter 'svm.path' must not be missing")
    }
  }
  
  if( any(classifyMethod==c("libsvm","svmlight","knn","nnet","ctreelibsvm")) ){
    classParameter <- switch(classifyMethod, 
                      "libsvm" = paste(paste("svm.kernel:",svm.kernel), 
                                   paste("svm.scale:",svm.scale), sep="; "), 

                      "svmlight" = paste(paste("svm.path:",svm.path), 
                                    paste("svm.options:",svm.options), sep="; "),          
                      "knn" = paste("knn.k:",knn.k),        
                      "nnet" = paste(paste("nnet.size:",nnet.size), 
                                paste("nnet.rang:",nnet.rang), 
                                paste("nnet.decay:",nnet.decay), 
                                paste("nnet.maxit:",nnet.maxit),
                                sep="; "), 
                      "ctreelibsvm" = paste(paste("svm.kernel:",svm.kernel), 
                                       paste("svm.scale:",svm.scale), sep="; "), 
    )    
  }else{
    classParameter <- "NULL"
  }
  
  addFeatureSet <- function(tmpData, featureFunction, featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit){    
    tmpModel = classify(tmpData, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)
    tmpPerformance =  tmpModel$totalPerformance    
    tmpNumber = length(tmpModel$features[[1]])
    tmpModel = c(featureFunction, featureParameter, ncol(tmpData)-1,
      classifyMethod, classParameter, cv)
    names(tmpModel) = c("Feature_Function","Feature_Parameter","Feature_Number",      
      "Model","Model_Parameter", "Cross_Validataion")    
    list(data=tmpData,model=tmpModel,performance=tmpPerformance) 
  }
    
  testFeatureSet = list()  
  
  if( any(featureMethod=="Binary") ){     
    tmpData = data.frame(featureBinary(seq,elements(ele.type)), classLable)  
    featureParameter <- paste("class:", ele.type)
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featureBinary", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)  
    if( ele.type=="aminoacid" ){      
      for(aa.type in group){
        tmpData = data.frame(featureBinary(seq,aaClass(aa.type)), classLable)  
        featureParameter <- paste("class:", aa.type)        
        testFeatureSet[[length(testFeatureSet)+1]] =
          addFeatureSet(tmpData, "featureBinary", featureParameter, classifyMethod=classifyMethod, cv=cv, 
          svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
          knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
          nnet.maxit=nnet.maxit)         
      }
    }
  }
  
  if( any(featureMethod=="CTD") ){
    tmpData = data.frame(featureCTD(seq,elements(ele.type)), classLable)  
    featureParameter <- paste("class:", ele.type)
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featureCTD", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)  
    if( ele.type=="aminoacid" ){      
      for(aa.type in group){
        tmpData = data.frame(featureCTD(seq,aaClass(aa.type)), classLable)  
        featureParameter <- paste("class:", aa.type)        
        testFeatureSet[[length(testFeatureSet)+1]] =
          addFeatureSet(tmpData, "featureCTD", featureParameter, classifyMethod=classifyMethod, cv=cv, 
          svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
          knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
          nnet.maxit=nnet.maxit)         
      }
    }
  }  
  
  if( any(featureMethod=="FragmentComposition") ){
    if ( missing(k) ) {
      stop("Parameter 'k' is needed for 'FragmentComposition'")
    }  
    tmpData = data.frame(featureFragmentComposition(seq,k,elements(ele.type)), classLable)  
    featureParameter <- paste(paste("class:", ele.type), paste("k: ", k), sep="; ")
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featureFragmentComposition", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)  
    if( ele.type=="aminoacid" ){      
      for(aa.type in group){
        tmpData = data.frame(featureFragmentComposition(seq,k,aaClass(aa.type)), classLable)  
        featureParameter <- paste(paste("class:", aa.type), paste("k:", k), sep="; ")       
        testFeatureSet[[length(testFeatureSet)+1]] =
          addFeatureSet(tmpData, "featureFragmentComposition", featureParameter, classifyMethod=classifyMethod, cv=cv, 
          svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
          knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
          nnet.maxit=nnet.maxit)         
      }
    }  
  }
    
  if( any(featureMethod=="GapPairComposition") ){
    if ( missing(g) ) { 
      stop("Parameter 'g' is needed for 'GapPairComposition'")
    }  
    tmpData = data.frame(featureGapPairComposition(seq,g,elements(ele.type)), classLable)  
    featureParameter <- paste(paste("class:", ele.type), paste("g:", g), sep="; ")
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featureGapPairComposition", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)  
    if( ele.type=="aminoacid" ){      
      for(aa.type in group){
        tmpData = data.frame(featureGapPairComposition(seq,g,aaClass(aa.type)), classLable)  
        featureParameter <- paste(paste("class:", aa.type), paste("g:", g), sep="; ")       
        testFeatureSet[[length(testFeatureSet)+1]] =
          addFeatureSet(tmpData, "featureGapPairComposition", featureParameter, classifyMethod=classifyMethod, cv=cv, 
          svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
          knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
          nnet.maxit=nnet.maxit)         
      }
    }  
  }    
    
  if( any(featureMethod=="CKSAAP") ){
    if ( missing(g) ) { 
      stop("Parameter 'g' is needed for 'CKSAAP'")
    }  
    tmpData = data.frame(featureCKSAAP(seq,g,elements(ele.type)), classLable)  
    featureParameter <- paste(paste("class:", ele.type), paste("g:", g), sep="; ")
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featureCKSAAP", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)  
    if( ele.type=="aminoacid" ){      
      for(aa.type in group){
        tmpData = data.frame(featureCKSAAP(seq,g,aaClass(aa.type)), classLable)  
        featureParameter <- paste(paste("class:", aa.type), paste("g: ", g), sep="; ")       
        testFeatureSet[[length(testFeatureSet)+1]] =
          addFeatureSet(tmpData, "featureCKSAAP", featureParameter, classifyMethod=classifyMethod, cv=cv, 
          svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
          knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
          nnet.maxit=nnet.maxit)         
      }
    }  
  }     
    
  if( any(featureMethod=="Hydro") & ele.type=="aminoacid" ){
    if ( missing(hydro.methods) ) { 
      stop("Parameter 'hydro.methods' is needed for 'Hydro'")
    }  
    for ( hydro.method in hydro.methods ){
      tmpData = data.frame(featureHydro(seq,hydro.method), classLable)  
      featureParameter <- paste("hydro.method:", hydro.method)
      testFeatureSet[[length(testFeatureSet)+1]] = 
        addFeatureSet(tmpData, "featureHydro", featureParameter, classifyMethod=classifyMethod, cv=cv, 
        svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
        knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
        nnet.maxit=nnet.maxit)
    }
  } 
  
  if( any(featureMethod=="ACH") & ele.type=="aminoacid" ){
    if ( missing(hydro.indexs) ) { 
      stop("Parameter 'hydro.indexs' is needed for 'ACH'")
    }  
    for ( hydro.index in hydro.indexs ){
      tmpData = data.frame(featureACH(seq,hydro.index), classLable)  
      featureParameter <- paste("hydro.index:", hydro.index)
      testFeatureSet[[length(testFeatureSet)+1]] = 
        addFeatureSet(tmpData, "featureACH", featureParameter, classifyMethod=classifyMethod, cv=cv, 
        svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
        knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
        nnet.maxit=nnet.maxit)
    }
  } 
  
  data(aa.index)
  index = sapply(aa.index, function(x){x$I})
  index = index[,apply(index,2,function(x){sum(is.na(x))==0})]
  indexName = colnames(index)
   
  if( any(featureMethod=="AAindex") & ele.type=="aminoacid" ){
    for ( aaindex.name in indexName ){
      tmpData = data.frame(featureAAindex(seq,aaindex.name), classLable)  
      featureParameter <- paste("aaindex.name:", aaindex.name)
      testFeatureSet[[length(testFeatureSet)+1]] = 
        addFeatureSet(tmpData, "featureAAindex", featureParameter, classifyMethod=classifyMethod, cv=cv, 
        svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
        knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
        nnet.maxit=nnet.maxit)
    }
  } 
  
  if( any(featureMethod=="ACI") & ele.type=="aminoacid" ){
    for ( aaindex.name in indexName ){
      tmpData = data.frame(featureACI(seq,aaindex.name), classLable)  
      featureParameter <- paste("aaindex.name:", aaindex.name)
      testFeatureSet[[length(testFeatureSet)+1]] = 
        addFeatureSet(tmpData, "featureACI", featureParameter, classifyMethod=classifyMethod, cv=cv, 
        svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
        knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
        nnet.maxit=nnet.maxit)
    }
  }   
  
  if( any(featureMethod=="ACF") & ele.type=="aminoacid" ){
    if ( missing(n) ) { 
      stop("Parameter 'n' is needed for 'ACF'")
    }  
    for ( aaindex.name in indexName ){
      tmpData = data.frame(featureACF(seq,n,aaindex.name), classLable)  
      featureParameter <- paste(paste("aaindex.name:", aaindex.name), paste("n:", n), sep="; ")      
      testFeatureSet[[length(testFeatureSet)+1]] = 
        addFeatureSet(tmpData, "featureACF", featureParameter, classifyMethod=classifyMethod, cv=cv, 
        svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
        knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
        nnet.maxit=nnet.maxit)
    }
  }   
  
  if( any(featureMethod=="PseudoAAComp") & ele.type=="aminoacid" ){
    if ( missing(d) | missing(w) ) { 
      stop("Parameter 'd' and 'w' is needed for 'PseudoAAComp'")
    }      
    tmpData = data.frame(featurePseudoAAComp(seq,d,w), classLable)  
    featureParameter <- paste(paste("d:", d), paste("w:", w), sep="; ")      
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featurePseudoAAComp", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)    
  }
  
  if( any(featureMethod=="PSSM") & ele.type=="aminoacid" ){
    if ( any(missing(start.pos), missing(stop.pos), missing(psiblast.path), missing(database.path)) ) { 
      stop("Parameter 'start.pos', 'stop.pos', 'psiblast.path' and 'database.path' is needed for 'PseudoAAComp'")
    }      
    tmpData = data.frame(featurePSSM(seq, start.pos, stop.pos, psiblast.path, psiblast.path), classLable)  
    featureParameter <- paste(paste("psiblast.path:", psiblast.path), paste("psiblast.path:", psiblast.path), sep="; ")      
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featurePSSM", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)    
  }  
  
  if( any(featureMethod=="DOMAIN") & ele.type=="aminoacid" ){
    if ( any(missing(hmmpfam.path), missing(pfam.path), missing(Evalue)) ) { 
      stop("Parameter 'hmmpfam.path', 'pfam.path' and 'Evalue' is needed for 'DOMAIN'")
    }
    domain = predictPFAM(seq, hmmpfam.path, pfam.path, Evalue)  
    tmpData = data.frame(featureDOMAIN(domain), classLable)  
    featureParameter <- paste(paste("hmmpfam.path:", hmmpfam.path), paste("pfam.path:", pfam.path), paste("Evalue:", Evalue), sep="; ")      
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featureDOMAIN", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)    
  }
  
  if( any(featureMethod=="BDNAVIDEO") & any(ele.type==c("rnaBase","dnaBase")) ){    
    tmpData = data.frame(featureBDNAVIDEO(seq), classLable)  
    featureParameter <- "NULL"      
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featureBDNAVIDEO", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)    
  }
  
  if( any(featureMethod=="DIPRODB") & any(ele.type==c("rnaBase","dnaBase")) ){    
    if ( any(missing(na.type), missing(na.strand), missing(diprodb.method), missing(diprodb.type)) ) { 
      stop("Parameter 'na.type', 'na.strand', 'diprodb.method' and 'diprodb.type' is needed for 'DOMAIN'")
    }  
    tmpData = data.frame(featureDIPRODB(seq, na.type, na.strand, diprodb.method, diprodb.type), classLable)  
    featureParameter <- paste(paste("na.type:", na.type), paste("na.strand:", na.strand), paste("diprodb.method:", diprodb.method), paste("diprodb.type:",diprodb.type), sep="; ")      
    testFeatureSet[[length(testFeatureSet)+1]] = 
      addFeatureSet(tmpData, "featureDIPRODB", featureParameter, classifyMethod=classifyMethod, cv=cv, 
      svm.kernel=svm.kernel, svm.scale=svm.scale, svm.path=svm.path, svm.options=svm.options, 
      knn.k=knn.k, nnet.size=nnet.size, nnet.rang=nnet.rang, nnet.decay=nnet.decay, 
      nnet.maxit=nnet.maxit)    
  }  
 
  summary <- t(sapply(testFeatureSet,function(x){
    c(x$model, format(x$performance,digits=2) )
  }))
  ord <- order(as.numeric(summary[,"acc"]), decreasing=T)
  summary <- summary[ord,]
  write.table(summary,fileName,quote=F,sep="\t",row.names=F)  
  testFeatureSet <- testFeatureSet[ord]
  testFeatureSet
}

