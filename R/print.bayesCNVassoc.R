print.bayesCNVassoc<-function(x, sig.CNVs=c(0.05,0.95), ...)
 {

  if (!inherits(x, "bayesCNVassoc")) 
     stop("object must be of class 'bayesCNVassoc'")

  N.groups <- x$N.groups
  v.median.ind.2 <- x$res.summary$v.median.ind.2
  v.median.select<-list()
  for (i in 1:N.groups)
   {
     v.median.select[[i]]<-v.median.ind.2[[i]][(v.median.ind.2[[i]][,4]>sig.CNVs[2] |    v.median.ind.2[[i]][,4]<sig.CNVs[1]),]
   }
  names(v.median.select)<-x$names.groups

  print(v.median.select)
}


