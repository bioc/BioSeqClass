hr <- function(seq, method, identity, cdhit.path){
    ## Filter homolog Reduction.
    ##
    ## seq - a list with one element for each protein/gene sequence. The elements 
    ##       are in two parts, one the description and the second a character 
    ##       string of the biological sequence.   
    ## identity - a numeric value ranged from 0 to 1. It is used as a 
    ##             maximum identity cutoff among input sequences.         
    ## method - a string for the method of homolog redunction.
    ##          "cdhit": http://www.bioinformatics.org/download.php/cd-hit/cd-hit-2007-0131.tar.gz
    ##                   it is suitable for filter full-length protein or gene sequences.
    ##          "aligndis": it is suitable for filter aligned seuqnces with same length.
    ## cdhit.path - a string for the path of cdhit program. eg: "/people/hongli/cd-hit/".
    ##
    ## Copyright 2008, Hong Li, all rights reserved.
    ##
    if( any(missing(seq), missing(method),missing(identity)) ){
      stop("Parameters seq, method and identity can not be missing or NULL")
    }  
    if( identity<0 | identity>1 ){
    	stop(paste("identity",identity,"is not supported. Has to be a numeric 
    	value ranged from 0 to 1"))
    }
    reducSeq = switch(method,
                       "cdhit" = cdhitHR(seq,identity,cdhit.path),
                       "aligndis" = aligndisHR(seq,identity),
                       stop(paste("method",method,"is not supported. Has to be cdhit or aligndis"))
                )
   return(reducSeq)           
}

## There is a litter bug in "psi-cd-hit-local.pl" in the cd-hit program package.   
## Please change the 
## "  $cmd = `$psi_blast -d ./$tmp_db $bl_para -i $seq_dir/$id -R $prof_dir/$id.prof -o $bl_dir/$id`;
##    $cmd = `$blastp -d ./$tmp_db $bl_para -i $seq_dir/$id -o $bl_dir/$id`; 
## "
## to
## "  $cmd = `$psi_blast -d $tmp_db $bl_para -i $seq_dir/$id -R $prof_dir/$id.prof -o $bl_dir/$id`;
##    $cmd = `$blastp -d $tmp_db $bl_para -i $seq_dir/$id -o $bl_dir/$id`; 
## "
cdhitHR <- function(seq, identity=0.3, cdhit.path){
    if( !file.exists(file.path(cdhit.path,"cd-hit")) | 
        !file.exists(file.path(cdhit.path,"psi-cd-hit.pl")) ){
      stop(paste("cdhit.path",cdhit.path,"is not corrected"))
    }
    names(seq) = sapply(seq,function(x){x$desc})
    perlName <- paste(tempfile("tempPerl"), "pl", sep=".")
    .pathPerl(perlName,.Platform$OS.type)    
    inFasta = tempfile("tempFasta") ;
    outFile = tempfile("tempOut") ;
    writeFASTA(seq,inFasta)
    write(paste("$inFasta='",inFasta,"';",sep=""), file = perlName, append = TRUE)     
    write(paste("$outFile='",outFile,"';",sep=""), file = perlName, append = TRUE)     
    write(paste("$identity=",identity,";"), file = perlName, append = TRUE) 
    write(paste("$cdhit='",cdhit.path,"';",sep=""), file = perlName, append = TRUE) 
    write(readLines( file.path(.path.package("BioSeqClass"),"scripts",
      "homologReduction_cdhit") ), file = perlName, append = TRUE) 
    .callPerl(perlName, .Platform$OS.type) 
    reducSeq = seq[as.vector(as.matrix(read.csv(outFile,header=F,sep="\t")))]     
    
    unlink(perlName)
    unlink(file.path(dirname(inFasta),dir(path=dirname(inFasta),pattern=basename(inFasta))),recursive=T) 
    return(reducSeq)
}


aligndisHR <- function(seq, identity=0.6){
## Ref: "Daniel Schwartz, Michael F. Chou, George M. Church. Predicting protein 
##      post-translational modifications using meta-analysis of proteome-scale 
##¡¡¡¡¡¡data sets. MCP."  
  D <- ceiling(length(unlist(strsplit(seq[[1]][["seq"]],split="")))*(1-identity))  
  tmpSeq2 <- sapply(seq,function(x){x[["seq"]]})
  names(tmpSeq2) <- sapply(seq,function(x){x[["desc"]]}) 
  for( j in 0:D ){
    for(i in names(tmpSeq2)){
      if( !is.na(tmpSeq2[i]) ){      
        tmp <- sapply(tmpSeq2[setdiff(names(tmpSeq2),i)], 
               function(x){distance(tmpSeq2[i],x)})
        tmpSeq2 <- tmpSeq2[setdiff(names(tmpSeq2),names(tmp)[tmp==j])]
      }
    }    
  }
  
  reducSeq <- list()                  
  for (x in 1:length(tmpSeq2) ){reducSeq[[x]]=list("desc"=names(tmpSeq2)[x],"seq"=tmpSeq2[x]) }  
  return(reducSeq)
}

distance <- function(seq1,seq2)
## Caclate the number of different residues between seq1 and seq2.
##
## seq1 - a string for the protein or gene sequence
## seq2 - a string for the protein or gene sequence.
##
{
  if( missing(seq1) | missing(seq2) ){
    stop("Parameters seq1 and seq2 can not be missing or NULL")
  }
  if( length(seq1)>1 | length(seq2)>1 ){
    stop(paste("Length of",seq1,"or",seq2,"must be 1"))
  }
  tmp1 <- unlist(strsplit(seq1,split=""))
  tmp2 <- unlist(strsplit(seq2,split=""))
  if( length(tmp1)!=length(tmp2) ){
    stop(paste(seq1,"and",seq2,"must have the same length"))
  }
  sum(tmp1!=tmp2)
}


getNegSite <- function(posSite, seq, aa){
  protein = unique(sapply(posSite,function(x){unlist(strsplit(x,split=":"))[1]}))
  negSite = unlist(sapply(protein,function(x){
              tmp = unlist(strsplit(seq[[x]][["seq"]],split=""))
              paste(x,(1:length(tmp))[tmp==aa],sep=":")
            }))
  negSite = setdiff(negSite,posSite)          
  negSite
}

getTrain <- function(seqfile, posfile, aa, w, identity, balance = T){
    seq = readFASTA(seqfile, strip.desc=TRUE)
    names(seq) = sapply(seq,function(x){x[["desc"]]})    
    tmp = as.matrix(read.csv(posfile, sep = "\t",header=F))
    tmp = tmp[apply(tmp, 1, function(x) {
        as.numeric(x[2]) - w > 0 & as.numeric(x[2]) + w <= length(unlist(strsplit(seq[[x[1]]][["seq"]], 
            split = "")))
    }), ]
    posSeq = list()
    for (x in 1:nrow(tmp) ){posSeq[[x]]=list("desc"=paste(tmp[x, 1], as.numeric(tmp[x, 2]), sep = ":"),
                                             "seq"=substr(seq[[tmp[x, 1]]][["seq"]], as.numeric(tmp[x, 2]) - w, as.numeric(tmp[x, 2]) + w) ) }        
    print(paste("Positive Site:", length(posSeq)))
    print(paste("Positive Protein:", length(unique(sapply(sapply(posSeq,function(x){x[["desc"]]}), 
        function(x) {
            unlist(strsplit(x, split = ":"))[1]
        })))))
    posSeq = hr(posSeq, method = "aligndis", identity = identity)
    print(paste("Positive Site After Homolog Reduction:", length(posSeq)))
    protein = unique(sapply(sapply(posSeq,function(x){x[["desc"]]}), function(x) {
        unlist(strsplit(x, split = ":"))[1]
    }))
    print(paste("Positive Protein After Homolog Reduction:", 
        length(protein)))
    negSite = getNegSite(sapply(posSeq,function(x){x[["desc"]]}), seq, aa)
    tmp = t(sapply(negSite, function(x) {
        unlist(strsplit(x, split = ":"))
    }))
    tmp = tmp[apply(tmp, 1, function(x) {
        as.numeric(x[2]) - w > 0 & as.numeric(x[2]) + w <= length(unlist(strsplit(seq[[x[1]]][["seq"]], 
            split = "")))
    }), ]
    negSeq = list()
    for (x in 1:nrow(tmp) ){negSeq[[x]]=list("desc"=paste(tmp[x, 1], as.numeric(tmp[x, 2]), sep = ":"),
                                             "seq"=substr(seq[[tmp[x, 1]]][["seq"]], as.numeric(tmp[x, 2]) - w, as.numeric(tmp[x, 2]) + w) ) }        

    tmp1=sapply(posSeq,function(x){x[["seq"]]})
    tmp2=sapply(negSeq,function(x){x[["seq"]]})    
    negSeq = negSeq[sapply(tmp2,function(x){  sum(x==tmp1)==0 })]
    print(paste("Negative Site:", length(negSeq)))
    print(paste("Negative Protein:", length(unique(sapply(sapply(negSeq,function(x){x[["desc"]]}), 
        function(x) {
            unlist(strsplit(x, split = ":"))[1]
        })))))
    if (balance) {
        if( length(posSeq)*2<length(negSeq) ){
          negSeq = sample(negSeq, length(posSeq)*2 )          
        }
        negSeq = hr(negSeq, method = "aligndis", identity = identity)
        if( length(posSeq)<length(negSeq) ){
          negSeq = sample(negSeq, length(posSeq))
        }else{
          print( "The number of Negative site is smaller than the Positive site")
        }
        print(paste("Negative Site After Homolog Reduction and Random Selection:", 
            length(negSeq)))
        print(paste("Negative Protein After Homolog Reduction and Random Selection:", 
            length(unique(sapply(sapply(negSeq,function(x){x[["desc"]]}), function(x) {
                unlist(strsplit(x, split = ":"))[1]
            })))))
    }
    else {
        negSeq = hr(negSeq, method = "aligndis", identity = identity)
        print(paste("Negative Site After Homolog Reduction:", 
            length(negSeq)))
        print(paste("Negative Protein After Homolog Reduction:", 
            length(unique(sapply(sapply(negSeq,function(x){x[["desc"]]}), function(x) {
                unlist(strsplit(x, split = ":"))[1]
            })))))
    }
    list(pos = posSeq, neg = negSeq)
}
