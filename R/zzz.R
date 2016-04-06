
# from http://stackoverflow.com/a/34031214/470769
Sys.which2 <- function(cmd) {
    stopifnot(length(cmd) == 1)
    if (.Platform$OS.type == "windows") {
        suppressWarnings({
            pathname <- shell(sprintf("where %s 2> NUL", cmd), intern=TRUE)[1]
        })
        if (!is.na(pathname)) return(stats::setNames(pathname, cmd))
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
