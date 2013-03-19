\name{featureSSC}
\alias{featureSSC}
\alias{getDSSP}
\alias{predictPROTEUS}
\title{Feature Coding by secondary structure}
\description{
  It is suitable for peptides with odd residues and the 
  central residue has important role.
}
\usage{
  featureSSC(secondaryStructure, confidenceScore)

  # secondary structure from DSSP database
  getDSSP(pdb)
  # Protein secondary structure prediction
  predictPROTEUS(seq,proteus2.organism="euk")
}
\arguments{      
  \item{secondaryStructure}{a string vector for the protein secondary structure.
    It is consisted of three kinds of secondary structures: H = Helix, 
    E = Beta Strand, C = Coil.}
  \item{confidenceScore}{a string vector for the confidence score of secondary 
    structure prediction (0-9, 0 = low, 9 = high).}      
  \item{pdb}{a string vector for the name of pdb structure. (e.g. "43ca")} 
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.}     
  \item{proteus2.organism}{a string for the organism of proteus2 program. This 
    must be one of the strings "gram-", "gram+", "euk".}      
}
\details{  
\code{\link{featureSSC}} codes for the secondary structure of the central 
  residue of peptides. It is suitable for peptides with odd residues and the 
  central residue has important role.
  
  \code{\link{getDSSP}} returns a vector of secondary structure extracted from 
  DSSP database (\url{http://swift.cmbi.ru.nl/gv/dssp/}).
  
  \code{\link{predictPROTEUS}} predicts secondary structure based on protein 
  sequence using following methods : "PROTEUS2", "PSIPRED", "JNET", "TRANSSEC", 
  "JURY-OF-EXPERTS PREDICTION".  Parameter "proteus2.organism" can be "gram-" 
  for "Gram negative prokaryote", "gram+" for "Gram positive prokaryote", "euk" 
  for "Eukaryote". It returns.....
}
\author{Hong Li}
\examples{
if(interactive()){  
  file = file.path(path.package("BioSeqClass"), "example", "acetylation_K.fasta")  
  tmp = readAAStringSet(file) 
  proteinSeq = as.character(tmp)
  
  DSSP1 = getDSSP(c("108l","43ca"))
  DSSP2 = getDSSP(c("108l","43ca","aaaa")) 
  
  ## Predict protein secordary strucutre   
  PROTEUS = predictPROTEUS(proteinSeq[1:2],proteus2.organism="euk") 
   
  ## Use general feature conding functions to codes protein secordary strucutre 
  secondaryStructure = sapply(PROTEUS,function(x){paste(x[["PROTEUS2"]]$SecondaryStructure,collapse="")})
  confidenceScore = sapply(PROTEUS,function(x){paste(x[["PROTEUS2"]]$ConfidenceScore,collapse="")})  
  SSCTD = featureCTD(secondaryStructure, class=list("H"="H","E"="E","C"="C"))          
  # Codes for peptides which have equal length and their central residues are important
  secondaryStructure = sapply(PROTEUS,function(x){sub.seq(paste(x[["PROTEUS2"]]$SecondaryStructure,collapse=""), 1, 11)})
  confidenceScore = sapply(PROTEUS,function(x){sub.seq(paste(x[["PROTEUS2"]]$ConfidenceScore,collapse=""), 1, 11)})  
 
  SS1 = featureSSC(secondaryStructure, confidenceScore)  
}
}
