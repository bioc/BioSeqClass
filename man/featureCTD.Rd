\name{featureCTD}
\alias{featureCTD}
\title{Feature Coding by composition, transition and distribution}
\description{
  Sequences are coded based on their composition, transition and distribution.
}
\usage{  
  featureCTD(seq,class=elements("aminoacid"))    
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.} 
  \item{class}{a list for the class of biological properties. It can 
    be produced by \code{\link{elements}} and \code{\link{aaClass}}.}   
}
\details{     
  \code{\link{featureCTD}} returns a matrix with M+M*(M-1)/2+M*5 columns. Each row 
  represented features of one sequence coding by a M+M*(M-1)/2+M*5 dimension 
  numeric vector. Three kinds of coding: composition (C), transition (T) and
  distribution (D) are used. C is the number of amino acids of a particular 
  property (such as hydrophobicity) divided by the total number of amino acids. 
  T characterizes the percent frequency with which amino acids of a particular 
  property is followed by amino acids of a different property. D measures
  the chain length within which the first, 25, 50, 75 and 100% of the amino 
  acids of a particular property is located respectively.  
}
\author{Hong Li}
\examples{
if(interactive()){  
  file = file.path(path.package("BioSeqClass"), "example", "acetylation_K.fasta")  
  library(Biostrings)
  tmp = readAAStringSet(file)

  proteinSeq = as.character(tmp)
  
  CTD1 = featureCTD(proteinSeq, class=elements("aminoacid") )
  CTD2 = featureCTD(proteinSeq, class=aaClass("aaV") )
}
}
