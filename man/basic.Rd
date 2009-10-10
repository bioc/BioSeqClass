\name{basic}
\alias{basic}
\alias{elements}
\alias{aaClass}
\alias{pwm}
\alias{.pathPerl}
\alias{.callPerl}
\alias{sub.seq} 
\alias{dssp.ss}
\alias{aa.index}
\alias{PROPERTY}
\alias{DiProDB}
\title{Assistant Functions}
\description{
  Assistant functions including read/write files, invoke perl programs, and so on.
}
\usage{
  ## Elements and groups of base and amino acid
  elements(ele.type)
  aaClass(aa.type)    

  pwm(seq,class=elements("aminoacid"))
  
  sub.seq(seq, start, stop)  
  
  .pathPerl(perlName, os)
  .callPerl(perlName, os)  
  
  data(dssp.ss)
  data(aa.index)
  data(PROPERTY)
  data(DiProDB)   
}
\arguments{
  \item{ele.type}{a string for the type of biological sequence. This must be 
    one of the strings "rnaBase", "dnaBase", "aminoacid" or "aminoacid2".}   
  \item{aa.type}{a string for the group of amino acids. This must be one of the
    strings "aaH", "aaV", "aaZ", "aaP", "aaF", "aaS" or "aaE".}  
  \item{seq}{a string vector for the protein or gene sequences.} 
  \item{class}{a list for the class of biological properties. It can 
    be produced by \code{\link{elements}} and \code{\link{aaClass}}.} 
  \item{start}{a positive integer vector indicating the first element of subsequence.} 
  \item{stop}{a positive integer vector indicating the last element of subsequence.} 
  \item{perlName}{a character string for the name of perl program.} 
  \item{os}{character string, giving the Operating System (family) of the computer.}   
}
\details{                        
  \code{\link{pwm}} returns a M*N position weight matrix (PWM) of input sequences. 
  M is the number of elements given by parameter "class". N is the length of each 
  sequence. Each row is a kind of element, and each column is a position. The 
  input sequences must have equal length.

  \code{\link{.pathPerl}} write the path of Perl to perl program file.
  
  \code{\link{.callPerl}} invoke Perl program via R.
  
  \code{\link{sub.seq}} extract subsequences by given start and stop postions.
    
  \code{\link{dssp.ss}} is a vector storing the secondary structure data 
    from DSSP database (\url{http://swift.cmbi.ru.nl/gv/dssp/}).
    
  \code{\link{aa.index}} is a list storing the properties of amino acids 
    from AAIndex database (\url{http://www.genome.jp/aaindex}).
    
  \code{\link{PROPERTY}} is a list sotring the properties of dinucleotide 
    from B-DNA-VIDEO PROPERTY database (\url{http://wwwmgs.bionet.nsc.ru/mgs/systems/bdnavideo/}).  
    
  \code{\link{DiProDB}} is a list sotring the conformational and thermodynamic
   dinucleotide properties from DiProDB database (\url{http://diprodb.fli-leibniz.de/}).  
}
\author{Hong Li}
\examples{  
  ## use readFASTA()/writeFASTA() to read/write FASTA file.
  library(Biostrings)
  file = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.fasta")
  tmp = readFASTA(file)
  seq <- sapply(tmp,function(x){x[["seq"]]})
  names(seq) <- sapply(tmp,function(x){x[["desc"]]}) 
  
  ## get the sequence fragment  
  sub1 = sub.seq(seq,rep(1,length(seq)),rep(4,length(seq)))
  tmp = round(sapply(seq,function(x){runif(1,min=1,
    max=length(unlist(strsplit(x,split="")))-5)}))
  sub2 = sub.seq(seq,tmp,tmp+5)
  
  ## cacluate the position weight matrix for sequences with equal length
  m = BioSeqClass:::pwm(sub1)
  
  ## load data: dssp.ss
  data(dssp.ss)
  ## see the data in dssp.ss
  dssp.ss[1:5]
}

