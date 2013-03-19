\name{featureHydro}
\alias{featureHydro}
\alias{featureACH}
\title{Feature Coding by hydrophobicity}
\description{
  Protein sequences are coded based on their hydrophobicity.
}
\usage{
  featureHydro(seq,hydro.method="SARAH1")
  featureACH(seq,hydro.index="hydroE")
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.} 
  \item{hydro.method}{a string for the method of coding protein hydrophobic effect. 
    This must be one of the strings "kpm" or "SARAH1".}
  \item{hydro.index}{a string for the method of coding protein hydrophobic effect. 
    This must be one of the strings "hydroE", "hydroF" or "hydroC".}
}
\details{   
  \code{\link{featureHydro}} returns a matrix measuring the hydrophobic effect. 
  Parameter "hydro.method" supported following coding methods:
  "kpm": use a numeral to indicating the hydrophobic effect of amino acid. Each 
         sequence is coded by a N dimension numeric vector.
  "SARAH1": use a 5 dimension 0-1 vector to indicating the hydrophobic effect of 
            amino acid. Each sequence is coded by a 5*N dimension 0-1 vector.
  
  \code{\link{featureACH}} returns a matrix with (N-1)/2 columns. N is the 
  length of input sequence, andis N must be odd. Central residue of all 
  windows are the central residue of input sequence. Each column is
  the average cumulative hydrophobicity over a sliding window.  
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(path.package("BioSeqClass"), "example", "acetylation_K.pos40.pep")
  seq = as.matrix(read.csv(file,header=F,sep="\t",row.names=1))[,1]
   
  H1 = featureHydro(seq,"kpm")
  H2 = featureHydro(seq,"SARAH1")

  H3 = featureACH(seq,hydro.index="hydroE")
  H3 = featureACH(seq,hydro.index="hydroF")
  H3 = featureACH(seq,hydro.index="hydroC")
}
}
