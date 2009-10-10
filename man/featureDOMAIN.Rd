\name{featureDOMAIN}
\alias{featureDOMAIN}
\alias{predictPFAM}
\title{Feature Coding by doamin organization}
\description{
  Protein sequences are coded based on their domains.
}
\usage{
  featureDOMAIN(domain)
  
  # Protein Pfam domain prediction
  predictPFAM(seq, hmmpfam.path, pfam.path, Evalue=10^-5)  
}
\arguments{    
 \item{domain}{a list of protein domains. It can be produced by function 
    \code{\link{predictPFAM}}.}    
  \item{seq}{a string vector for the protein, DNA, or RNA sequences.}   
  \item{hmmpfam.path}{a string for the path of hammpfam program in HMMER. 
    hammpfam will be employed to predict domains using models in Pfam database.}
  \item{pfam.path}{a string for the path of pfam domain database.}
  \item{Evalue}{a numeric value for the E-value cutoff of perdicted Pfam domain.}  
}
\details{  
  \code{\link{featureDOMAIN}} uses Pfam domains to code 0-1 feature vector.
  
  \code{\link{predictPFAM}} predict Pfam domains by hmmpfam program. It returns 
  a list, each element is a vector which denotes the domain composition of a 
  protein.   
}
\author{Hong Li}
\examples{
if(interactive()){

}
}
