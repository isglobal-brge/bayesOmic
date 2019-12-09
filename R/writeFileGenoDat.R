
writeFileGenoDat<-function(gg, allesum=FALSE, caco)
{
  if (missing(caco))
   stop("caco argument is required")

  ggin<-as.vector(gg)

  # pending: write check for correct gg entries
  lev<-as.raw(c(0,1,2,3))
  nprobes<-NCOL(gg)
  nsub<-NROW(gg)
  nlev<-length(lev)
  
  if(allesum)
   {
    if(is.null(caco))
	  stop("caco variable has no entries \n")
	
    if(is.null(levels(caco)))
	  stop("caco variable must be factor \n")

    levcaco<-levels(caco)
    numlevcaco<-length(levcaco)
	levcacoInt<-1:numlevcaco
    cacoInt<-sapply(caco,function(x) levcacoInt[x==levcaco])
   }
  else 
   {
    # irrelevant integers
    cacoInt<-0
    levcaco<-0
    numlevcaco<-0
   }
  
  out<-.C("writeGenoDat", 
           as.raw(ggin), 
           as.integer(nprobes), 
           as.integer(nsub), 
           as.raw(lev), 
           as.integer(nlev), 
           as.integer(allesum), 
           as.integer(cacoInt), 
           as.integer(levcacoInt), 
           as.integer(numlevcaco), 
           outInv=as.double(rep(0,nprobes)), 
           outAleleSum=as.double(rep(0,numlevcaco*nprobes)), 
           outNoMissCount=as.double(rep(0,numlevcaco*nprobes)),
           PACKAGE="bayesGen")
          

  ans<-list(out$outInv, 
            t(matrix(out$outAleleSum, nrow=numlevcaco, ncol=nprobes)), 
            t(matrix(out$outNoMissCount, nrow=numlevcaco, ncol=nprobes))
           )
  ans
}
