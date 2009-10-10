\name{featureDIPRODB}
\alias{featureDIPRODB}
\title{Feature Coding by Dinucleotide Property}
\description{
  Sequences are coded by conformational or thermodynamic dinucleotide property from 
  DiProDB database.
}
\usage{
  featureDIPRODB(seq, na.type="all", na.strand="all", diprodb.method="all", 
    diprodb.type="all")
}
\arguments{    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.}   
  \item{na.type}{a string for nucleic acid type. It must be "DNA", "DNA/RNA", "RNA", 
    or "all".}
  \item{na.strand}{a string for strand information. It must be "double", "single", 
    or "all".}
  \item{diprodb.method}{a string for mode of property determination. It can be 
    "experimental", "calculated", or "all".}
  \item{diprodb.type}{a string for property type. It can be "physicochemical", 
    "conformational", "letter based", or "all".}
}
\details{  
  \code{\link{featureDIPRODB}} returns a matrix with 122 columns. Each column is
  the mean of conformational or thermodynamic dinucleotide property from 
  DiProDB database (\url{http://diprodb.fli-leibniz.de}).  
}
\author{Hong Li}
\examples{
if(interactive()){
  file = file.path(.path.package("BioSeqClass"), "example", "test.rna")    
  rna = as.matrix(read.csv(file,header=F,sep="\t"))[,1]
  
  DIPRODB1 = featureDIPRODB(rna)
  DIPRODB2 = featureDIPRODB(rna, na.type="RNA")
}
}
