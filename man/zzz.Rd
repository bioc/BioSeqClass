\name{zzz}
\alias{.onLoad}
\title{Load packages and data}
\description{
  This functions load depended R packages and imports default data into 
  global "options".
}
\usage{
  .onLoad(libname, pkgname)
}
\arguments{
  \item{libname}{a character string giving the library directory where the 
    package defining the namespace was found.}
  \item{pkgname}{a character string giving the name of the package.}  
}
\details{
  After loading, loadNamespace looks for a hook function named .onLoad and runs 
  it before sealing the namespace and processing exports.
}
\author{Hong Li}