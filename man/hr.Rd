\name{hr}
\alias{hr}
\alias{cdhitHR}
\alias{aligndisHR}
\alias{getTrain}
\alias{getNegSite}
\alias{distance}
\title{Homolog Reduction}
\description{
  Filter homolog sequences by sequence similarity.
}
\usage{
  hr(seq, method, identity, cdhit.path)
  
  cdhitHR(seq, identity=0.3, cdhit.path)
  aligndisHR(seq, identity=0.6)
  distance(seq1,seq2)
  
  getTrain(seqfile, posfile, aa, w, identity, balance=T)  
  getNegSite(posSite, seq, aa)
}
\arguments{
  \item{seq}{a list with one element for each protein/gene sequence. The 
    elements are in two parts, one the description ("desc") and the second is a 
    character string of the biological sequence ("seq").} 
  \item{identity}{a numeric value ranged from 0 to 1. It is used as a maximum 
    identity cutoff among input sequences. } 
  \item{method}{a string for the method of homolog redunction. This must be 
    one of the strings "cdhit" or "aligndis".} 
  \item{cdhit.path}{a string for the path of cdhit program directory. eg: 
    "/people/hongli/cd-hit". It is necessary when method="cdhit".} 
  \item{seq1}{a string for the protein or gene sequence.}    
  \item{seq2}{a string for the protein or gene sequence. seq1 and seq2 must 
    have same length.}
  \item{seqfile}{a string for the name of FASTA file.} 
  \item{posfile}{a string for the name of file which contains the positive site 
    dataset. It has two columns: 1st column is the protein name; 2st column is 
    the positive site. Protein name should be consistent with the name used in 
    seqfile.} 
  \item{aa}{a character for the interested amino acid. eg: "C".}
  \item{w}{an integer for the window size of flanking peptide sequence. Window 
    size is 2*w+1, and the central residues are the positive sites in posfile.}
  \item{balance}{a logical value indicating whether negative sites will be random
    selected to have the same number with positive sites.}
  \item{posSite}{a string vector for the positive sites. It is consisted of protein
    description and positive site, eg: "P278168:952".} 
}
\details{  
  \code{\link{hr}} employs \code{\link{cdhitHR}} and \code{\link{aligndisHR}} to
  filter homolog sequences. It supported following methods:
  
  "cdhit": Use cd-hit program to quickly filter sequences by given identity. 
           It is designed to filter full-length protein or gene sequences.
           "formatdb" and "blastall" are required for running cd-hit program.         
           (http://www.bioinformatics.org/download.php/cd-hit/cd-hit-2007-0131.tar.gz
           or http://www.bioinformatics.org/download.php/cd-hit/cd-hit-2007-0131-win32.tar.gz)
  
  "aligndis": Use the number of different residues to meature the identity 
              between two sequences. 
              It is designed to filter aligned seuqnces with equal length.
  
  \code{\link{getTrain}} extract 2*w+1 flanking peptides of positive sites and 
    filter homolog sequences. Negative sites are non-positive sites in the same
    proteins. 
    
  \code{\link{distance}} calculate the number of positions with different residues
  between two sequences.
}
\value{
  \code{\link{hr}} return a list of reduced sequences.
}
\author{Hong Li}
\examples{
  distance("AABD","ACBD")
  distance("AABD","ECBD")
  if(interactive()){  
    file = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.fasta")
    library(Biostrings)
    seq = as.character(readAAStringSet(file))
    ## Homolog reduction of whole-length sequence by cd-hit
    # need cd-hit program;
    reducSeq50 = hr(seq, method="cdhit", identity=0.5, cdhit.path="/people/hongli/cd-hit")
    
    file = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.site")
    tmp = as.matrix(read.csv(file, sep="\t",header=F))
    logical = apply(tmp,1,function(x){ l=nchar(seq[x[1]]); (l>=as.numeric(x[2])+7 & as.numeric(x[2])-7>0) })
    fragment = sub.seq(seq[tmp[logical,1]], as.numeric(tmp[logical,2])-7, as.numeric(tmp[logical,2])+7)  
    ## Homolog reduction of short sequence fragment
    # It may be slow.
    reducSeq = hr(fragment, method="aligndis", identity=0.4)
    
    ## produce train set based on given positive sites and fasta sequences. 
    file = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.fasta")
    posfile = file.path(.path.package("BioSeqClass"), "example", "acetylation_K.site")
    ## "getTrain" integrate negative set construction and homolog reduction. It is designed for site level training data. 
    # It may be very slow.
    data = getTrain(file, posfile, aa="K", w=7, identity=0.4)
  }
}
