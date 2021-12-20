#' @keywords internal
.onAttach <- function(libname, pkgname){
  msg <- paste("Thank you for using our package. Currently we support using lm and glm. When using glm, please specify family (quasibinomial for binary outcome).")
  packageStartupMessage(msg)
  msg <- paste0('We politely ask you to cite our package using: citation("', pkgname, '"). It will be appreciated. ')
  packageStartupMessage(msg)
}