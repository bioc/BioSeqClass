\name{featureFragmentComposition}
\alias{featureFragmentComposition}
\title{Feature Coding by the composition of k-mer fragments}
\description{
  Sequences are coded based on the frequency of k-mer sequence fragments.
}
\usage{
  featureFragmentComposition(seq,k,class=elements("aminoacid"))
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.} 
  \item{k}{an integer indicating the length of sequence fragment (k>=1).} 
  \item{class}{a list for the class of biological properties. It can 
    be produced by \code{\link{elements}} and \code{\link{aaClass}}.}   
}
\details{     
  \code{\link{featureFragmentComposition}} returns a matrix with M\^k columns. 
  Each row represented features of one sequence coding by a M\^k dimension 
  numeric vector. Each column is the frequency of k-mer sequence fragment.
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  seq = as.matrix(read.csv(file,header=F,sep="\t",row.names=1))[,1]
  
  FC2 = featureFragmentComposition(seq,2,aaClass("aaS"))
  FC3 = featureFragmentComposition(seq,3,aaClass("aaS"))
}
}
