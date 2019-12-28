#' Get Parameters for a Bayesian model (SNP or CNV)
#' 
#' @aliases getParameters
#' @aliases print.parameters
#' @param x An object of class 'bayesSNPassoc' or 'bayesCNVassoc'
#' @method print parameters




getParameters <- function(x)
 {
  alpha<-x$res.summary$alpha
  loglambda<-x$res.summary$loglambda

  ans<-list(alpha=alpha, loglambda=loglambda)
  
  class(ans)<-"parameters"
  
  ans
 }

