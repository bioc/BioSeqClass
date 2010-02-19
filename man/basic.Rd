\name{basic}
\alias{basic}
\alias{elements}
\alias{aaClass}
\alias{pwm}
\alias{.pathPerl}
\alias{.callPerl}
\alias{sub.seq} 
\alias{dssp.ss}
\alias{aa.index}
\alias{PROPERTY}
\alias{DiProDB}
\title{Assistant Functions}
\description{
  Assistant functions including read/write files, invoke perl programs, and so on.
}
\usage{
  ## Elements and groups of base and amino acid
  elements(ele.type)
  aaClass(aa.type)    

  pwm(seq,class=elements("aminoacid"))
  
  .pathPerl(perlName, os)
  .callPerl(perlName, os)  
  
  data(dssp.ss)
  data(aa.index)
  data(PROPERTY)
  data(DiProDB)   
}
\arguments{
  \item{ele.type}{a string for the type of biological sequence. This must be 
    one of the strings "rnaBase", "dnaBase", "aminoacid" or "aminoacid2".}   
  \item{aa.type}{a string for the group of amino acids. This must be one of the
    strings "aaH", "aaV", "aaZ", "aaP", "aaF", "aaS" or "aaE".}  
  \item{seq}{a string vector for the protein or gene sequences.} 
  \item{class}{a list for the class of biological properties. It can 
    be produced by \code{\link{elements}} and \code{\link{aaClass}}.}   
  \item{perlName}{a character string for the name of perl program.} 
  \item{os}{character string, giving the Operating System (family) of the computer.}   
}
\details{     
  \code{\link{elements}} returns a list of basic elements of biological sequence. 
  Parameter "ele.type" supported following selection:  
  "rnaBase" - basic elements of RNA (ATCG).  
  "dnaBase" - basic elements of DNA (AUCG).
  "aminoacid" - 20 amino acides (RKEDQNWGASTPHYCVLIMF).
  "aminoacid2" - 20 amino acides and 1 pseudo amino acid "O" (RKEDQNWGASTPHYCVLIMFO). 
                 Unknown or uncomplete amino acides will be substituted by pseudo                  
                 amino acid.

  \code{\link{aaClass}} returns a list of amino acids groups depend on their 
  physical-chemical properties. Parameter "aa.type" supports following selection:
  "aaH" (hydrophobicity): Polar(RKEDQN), Neutral(GASTPHY), Hydrophobic(CVLIMFW)
  "aaV (normalized Van der Waals volume)": Small(GASCTPD), Medium(NVEQIL), Large(MHKFRYW)
  "aaZ" (polarizability): Low polarizability	(GASDT), Medium polarizability 
                          (CPNVEQIL), High polarizability(KMHFRYW)
  "aaP" (polarity): Low polarity	(LIFWCMVY), Neutral polarity (PATGS), High polarity	(HQRKNED)
  "aaF": Acidic	(DE), Basic (HKR), Polar (CGNQSTY), Nonpolar (AFILMPVW)
  "aaS": Acidic	(DE), Basic	(HKR), Aromatic (FWY), Amide (NQ), Small hydroxyl	(ST), 
         Sulfur-containing (CM), Aliphatic (AGPILV)
  "aaE": Acidic	(DE), Basic(HKR), Aromatic (FWY), Amide (NQ), Small hydroxyl (ST),
         Sulfur-containing (CM), Aliphatic 1 (AGP), Aliphatic 2	(ILV)
                   
  \code{\link{pwm}} returns a M*N position weight matrix (PWM) of input sequences. 
  M is the number of elements given by parameter "class". N is the length of each 
  sequence. Each row is a kind of element, and each column is a position. The 
  input sequences must have equal length.

  \code{\link{.pathPerl}} write the path of Perl to perl program file.
  
  \code{\link{.callPerl}} invoke Perl program via R.
    
  \code{\link{dssp.ss}} is a vector storing the secondary structure data 
    from DSSP database (\url{http://swift.cmbi.ru.nl/gv/dssp/}).
    
  \code{\link{aa.index}} is a list storing the properties of amino acids 
    from AAIndex database (\url{http://www.genome.jp/aaindex}).
    
  \code{\link{PROPERTY}} is a list sotring the properties of dinucleotide 
    from B-DNA-VIDEO PROPERTY database (\url{http://wwwmgs.bionet.nsc.ru/mgs/systems/bdnavideo/}).  
    
  \code{\link{DiProDB}} is a list sotring the conformational and thermodynamic
   dinucleotide properties from DiProDB database (\url{http://diprodb.fli-leibniz.de/}).  
}
\author{Hong Li}
\examples{  
  ## amino acids groups depend on their hydrophobicity
  aaClass("aaH")
  
  ## load data: dssp.ss
  data(dssp.ss)
  ## see the data in dssp.ss
  dssp.ss[1:5]
}

