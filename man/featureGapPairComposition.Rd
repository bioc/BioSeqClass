\name{featureGapPairComposition}
\alias{featureGapPairComposition}
\title{Feature Coding by g-spaced aminoacids/bases pairs}
\description{
  Sequences are coded based on the frequency of g-spaced aminoacids/bases pairs.
}
\usage{
  featureGapPairComposition(seq,g,class=elements("aminoacid"))  
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.} 
  \item{g}{an integer indicating the distance between two aminoacids/bases
    (g>=0).}  
  \item{class}{a list for the class of biological properties. It can 
    be produced by \code{\link{elements}} and \code{\link{aaClass}}.}   
}
\details{  
  \code{\link{featureGapPairComposition}} returns a matrix with M\^2 columns. 
  Each row represented features of one sequence coding by a M\^2 dimension 
  numeric vector. Each column is the frequency of g-spaced aminoacids/bases
  pair. featureFragmentComposition(seq,2) is same with featureGapPairComposition(seq,0). 
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  seq = as.matrix(read.csv(file,header=F,sep="\t",row.names=1))[,1]
  
  GPC0 = featureGapPairComposition(seq,0,elements("aminoacid"))
  GPC2 = featureGapPairComposition(seq,2,elements("aminoacid"))
}
}
