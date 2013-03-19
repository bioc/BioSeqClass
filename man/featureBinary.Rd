\name{featureBinary}
\alias{featureBinary}
\title{Feature Coding by Binary Vectors}
\description{
  Sequences are coded by binary vectors.
}
\usage{  
  featureBinary(seq,class=elements("aminoacid"))  
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.}
  \item{class}{a list for the class of biological properties. It can 
    be produced by \code{\link{elements}} and \code{\link{aaClass}}.}   
}
\details{    
  \code{\link{featureBinary}} returns a matrix with M*N columns. Each row 
  represented features of one sequence coding by a M*N dimension 0-1 vector. 
  Each base/amino acid is coded as a M dimension vetor. For example:  
  amino acid "A" is coded by "00000000000000000001"; base "T" is coded by 
  "0010". The input sequences must have equal length.  
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  seq = as.matrix(read.csv(file,header=F,sep="\t",row.names=1))[,1]
  
  BIN1 = featureBinary(seq,elements("aminoacid"))
  BIN2 = featureBinary(seq,aaClass("aaE"))
}
}
