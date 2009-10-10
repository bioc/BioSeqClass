\name{featurePseudoAAComp}
\alias{featurePseudoAAComp}
\title{Feature Coding by Pseudo Amino Acid Composiion}
\description{
  Protein sequences are coded by pseudo amino acid composiion.
}  
\usage{
  featurePseudoAAComp(seq,d,w=0.05)  
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.} 
  \item{d}{an integer used as paramter of \code{\link{featurePseudoAAComp}} 
    (d>=1). Coupling between amino acids X(i) and X(i+d) are considered as features.}
  \item{w}{a numeric value for the weight factor of sequence order effect in 
    \code{\link{featurePseudoAAComp}}. }
}
\details{     
  \code{\link{featurePseudoAAComp}} returns a matrix representing the pseudo 
  amino acid composiion. Each row represented features of one sequence coding 
  by a 20+d dimension numeric vector. The first 20 features indicates the 
  composition of 20 amino acids. The last d features indicates the coupling 
  between amino acids X(i) and X(i+d). Coupling value is cacluated by hydrophobicity, 
  hydrophilicity and mass of amino acids.  
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  seq = as.matrix(read.csv(file,header=F,sep="\t",row.names=1))[,1]
    
  PAC4 = featurePseudoAAComp(seq,4)    
}
}
