\name{featureAAindex}
\alias{featureAAindex}
\alias{featureACI}
\alias{featureACF}
\title{Feature Coding by physicochemical/biochemical properties in AAindex }
\description{
  Protein sequences are coded based on the physicochemical/biochemical properties 
  of amino acids in AAindex database.
}
\usage{
  featureAAindex(seq,aaindex.name="all")
  featureACI(seq,aaindex.name="all")
  featureACF(seq,n,aaindex.name="all") 
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.} 
  \item{aaindex.name}{a string for the name of physicochemical and biochemical 
    properties in AAindx.}
  \item{n}{an integer used as paramter of \code{\link{featureACF}} (1<=n<=L-2, 
    L is the the length of sequence). featureACF takes the auto-correlation 
    between fragment X(1)...X(L-m) and X(m+1)...X(L) (1<=m<=n) as features. }  
}
\details{  
  \code{\link{featureAAindex}} returns a matrix measuring the physicochemical 
  and biochemical properties of amino acids by AAindex (\url{http://www.genome.jp/aaindex}).
  If parameter aaindex.name="all", 
  all properties in AAindex will be considered, and each row represented the 
  features of one sequence coding by a 531*N dimension numeric vector. If parameter 
  aaindex.name is a name of property in AAindex, each row represented the features
  of one sequence coding by a N dimension numeric vector.
  
  \code{\link{featureACI}} returns a matrix with 531 columns, measuring the 
  average cumulative value of AAindex. N is the length of input sequence, and 
  N must be odd. Central residue of all windows are the central residue of input 
  sequence. Each column is the average cumulative AAindex over a sliding window.  
  
  \code{\link{featureACF}} returns a matrix with 531*n columns, measuring the 
  Auto-Correlation Function (ACF) of AAindex. If parameter aaindex.name is a name 
  of property in AAindex, each row represented  the features of one sequence coding 
  by a n dimension numeric vector.  
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  seq = as.matrix(read.csv(file,header=F,sep="\t",row.names=1))[,1]
 
  AI_all = featureAAindex(seq)
  AI_ANDN920101 = featureAAindex(seq,"ANDN920101")
  
  ACI_all = featureACI(seq)
  ACI_ANDN920101 = featureACI(seq,"ANDN920101")
  
  ACF_all_1 = featureACF(seq,1)
  ACF_ANDN920101_3 = featureACF(seq,3,"ANDN920101") 
}
}
