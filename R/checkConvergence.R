#' Check model convergence of coefficients from the Bayesian model
#' 
#' @aliases checkConvergence
#' @param x object of class 'bayesOmic'
#' @param parameter 
#' @param type
#' 
#' @export

checkConvergence <- function(x, parameter="alpha", type="density")
 {
  if (!inherits(x, "bayesOmic")) 
     stop("object must be of class 'bayesOmic'")

  type.plot <- charmatch(type, c("density", "traceplot", "running_mean"))
  if (is.na(type.plot))
    stop("'type' argument must be density, traceplot or running_mean")

  fun <- ifelse(type.plot==1, ggmcmc::ggs_density, 
                ifelse(type.plot==2, ggmcmc::ggs_traceplot, ggmcmc::ggs_running))
  
  param<-charmatch(parameter, c("alpha", "beta", "shared", "specific"))
 
  if (is.na(param))
   stop("'parameter' argument must be 'alpha', 'beta', 'shared' or 'specific'")

  nn<-colnames(x$res[[1]])

  if (param==1) {
    o<-grep("alpha", nn)
    S <- ggmcmc::ggs(x$res[,o])
    fun(S)
  }
  else if (param==2) {
    o<-grep("beta", nn)
    S <- ggmcmc::ggs(x$res[,o])
    fun(S)
  }
  else if (param==3){
    o<-grep("^u", nn)
    S <- ggmcmc::ggs(x$res[,o])
    fun(S)
  }
  else if (param==4){
    o<-grep("^v", nn)
    S <- ggmcmc::ggs(x$res[,o])
    fun(S)
  }
}

