getSpecific <- function(x){
  ss <- x$res.summary
  rownames(ss) <- rep(x$names.CNVs, 3)
  sig <- apply(x$res.summary,1, function(x) x[1]<0 & x[3]<0 | x[1]>0 & x[3]>0)
  n <- x$N.cnvs
  out <- list()
  for (i in 1:x$N.groups) {
    sel <- (n*(i-1) + 1):(n*i)
    out[[i]] <- ss[sel,][sig[sel],]  
  }
  names(out) <- x$names.groups
  out
}
  