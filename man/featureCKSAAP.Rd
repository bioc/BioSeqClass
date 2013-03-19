\name{featureCKSAAP}
\alias{featureCKSAAP}
\title{Feature Coding by k-spaced Aminoacids/Base Pairs}
\description{
  Protein sequences are coded based on the frequency of k-spaced aminoacids/base pairs.
}
\usage{
  featureCKSAAP(seq,g,class=elements("aminoacid"))  
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.} 
  \item{g}{an integer indicating the distance between two aminoacids/bases
    (g>=0).}  
  \item{class}{a list for the class of biological properties. It can 
    be produced by \code{\link{elements}} and \code{\link{aaClass}}.}   
}
\details{  
  \code{\link{featureCKSAAP}} returns a matrix with (g+1)*M\^2 columns. Each 
  row represented features of one sequence coding by a (g+1)*M\^2 dimension 
  numeric vector. Each column is the number of k-spaced aminoacids/base
  pair (0<=k<=g).  
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  seq = as.matrix(read.csv(file,header=F,sep="\t",row.names=1))[,1]
    
  CKSAAP0 = featureCKSAAP(seq,0,elements("aminoacid"))
  CKSAAP2 = featureCKSAAP(seq,2,elements("aminoacid"))  
}
}
