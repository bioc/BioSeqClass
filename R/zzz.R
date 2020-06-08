
# from http://stackoverflow.com/a/34031214/470769
Sys.which2 <- function(cmd) {
    stopifnot(length(cmd) == 1)
    if (.Platform$OS.type == "windows") {
        suppressWarnings({
            pathname <- shell(sprintf("where %s 2> NUL", cmd), intern=TRUE)[1]
        })
        if (!is.na(pathname)) return(dQuote(stats::setNames(pathname, cmd)))
    }
    Sys.which(cmd)
}

.javaExecutable <- function() Sys.which2("java")



.onLoad <- function(libname, pkgname) {
  # add vigette to menu; support Windows only
  if(.Platform$OS.type == "windows" && interactive()
      && .Platform$GUI == "Rgui"){
          addVigs2WinMenu("BioSeqClass")
  }
  data("dssp.ss", "aa.index", "PROPERTY", "DiProDB", package="BioSeqClass")
}
.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.13")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
