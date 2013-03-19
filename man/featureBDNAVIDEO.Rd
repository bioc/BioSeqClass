\name{featureBDNAVIDEO}
\alias{featureBDNAVIDEO}
\title{Feature Coding by DNA/RNA property}
\description{
  DNA/RNA Sequences are coded with DNA or RNA property from B-DNA-VIDEO database. 
}
\usage{
  featureBDNAVIDEO(seq)
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.} 
}
\details{  
  \code{\link{featureBDNAVIDEO}} returns a matrix with 38 columns. Each column is
  the mean of DNA or RNA property from B-DNA-VIDEO database 
  (\url{http://wwwmgs.bionet.nsc.ru/mgs/systems/bdnavideo/}).  
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(path.package("BioSeqClass"), "example", "test.rna")    
  rna = as.matrix(read.csv(file,header=F,sep="\t"))[,1]
  
  BDNAVIDEO = featureBDNAVIDEO(rna)
}
}
