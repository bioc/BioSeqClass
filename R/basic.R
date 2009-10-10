pwm <- function(seq,class=elements("aminoacid")){
## PWM: position weight matrix. 
  if( length(unique(sapply(seq,function(x){length(unlist(strsplit(x,split="")))})))>1 ){
    stop("Sequences in seq must have equal length")
  }
  class2 = rep(names(class),sapply(class,length) )
  names(class2) = unlist(class)
  aaM = sapply(seq,function(x){ unlist(strsplit(x,split="")) })
  apply(aaM,1,function(x){ 
               y=rep(0,length=length(class))
               names(y) = names(class)
               tmp=table(class2[x]) ; tmp = tmp/sum(tmp)                               
               y[names(tmp)]=tmp
               y
  })
}

.callPerl <- function(perlName, os){
    if(os == "unix"){
        system(paste("chmod +x", perlName))
        system(perlName)
    }else if(os == "windows"){
        perlName <- gsub("/", "\\\\", perlName)
        system(paste("perl", perlName))
    }else{
        stop(paste("Do not know who to run perl under ", os))
    }
}

.pathPerl <- function(perlName, os){
    if(os == "unix"){
        perlBin <- system("which perl", intern = TRUE)
        if(length(perlBin) > 1){
            stop("Perl is not available!")
        }
        write(paste("#!", perlBin, "\n\n", sep = ""), file = perlName)
    }else if(os == "windows"){
        write("#!/usr/bin/perl -w\n\n", file = perlName)
    }
}

sub.seq <- function(seq, start, stop){
  if( any(missing(seq), missing(start),missing(stop)) ){
    stop("Parameters seq, start and stop can not be missing or NULL")
  }
  if( any(length(seq)!=length(start), length(seq)!=length(stop), 
    length(stop)!=length(start)) ){
    stop("Parameters seq, start and stop must have equal length")
  }  
  if( sum(start<=0)>0  | sum(stop<=0)>0 ){
    stop("Parameters start and stop must be position integer")
  }
  if( is.null(names(seq)) ){
    names(seq) = paste("seq",1:length(seq),sep="")
  }
  subSeq = sapply(1:length(seq),function(x){ 
             tmp = length(unlist(strsplit(seq[x],split="")))
             if( stop[x]>tmp ){
               stop("Parameters stop can not be larger than the sequence length")
             }
             substr(seq[x],start[x],stop[x]) 
           })
  names(subSeq) = paste(names(seq),start,stop,sep=":")
  subSeq
}


