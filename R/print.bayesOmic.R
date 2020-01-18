#' print.bayesOmicAssoc
#' 
#' @param x object of class 'bayesOmic'
#' @param x 
#' @param ...
#' 
#' @S3method print bayesOmic

print.bayesOmic <- function(x, ...)
 {
  alpha<-x$res.summary$alpha
  beta<-x$res.summary$beta
  cat("\n Intercepts (alpha): \n")
  print(alpha)
  cat("\n Coefficients of shared components (beta): \n")
  print(beta)
  cat("\n Use 'getSpecific()' and 'getShared()' functions to get specific or shared components, respectively. \n")
 }

