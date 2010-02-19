.onLoad <- function(libname, pkgname) {
  # add vigette to menu; support Windows only
  if(.Platform$OS.type == "windows" && interactive()
      && .Platform$GUI == "Rgui"){
          addVigs2WinMenu("BioSeqClass")
  } 
  data("dssp.ss", "aa.index", "PROPERTY", "DiProDB", package="BioSeqClass")
}
