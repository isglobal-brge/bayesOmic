makeHeatmap <-
function(x, quantiles=c(.2,.5,.6,.8), ...)
 {

   if (!inherits(x, "bayesSNPassoc")) 
     stop("object must be of class 'bayesSNPassoc'")

   Nvar<-x$Nvar
   Ngroups<-x$Ngroups
   temp<-unlist(x$res.summary$pi.mean)
   pi.mean<-matrix(temp, ncol=Ngroups, nrow=Nvar)
   colnames(pi.mean)<-x$names.groups

   ncolor<-length(quantiles)+1
   heatmap(pi.mean, breaks=c(-1000, quantile(pi.mean, quantiles), 1000), col=brewer.pal(ncolor, "YlOrBr"), cexRow=.6, cexCol=.8, scale="none", margins=c(12,10), main=expression(paste("logit(",pi[ij],")",sep="")), ...)

   cc<-round(quantile(pi.mean,c(.2,.5,.6,.8)),4)

   my.leg<-rep(NA, ncolor)
   my.leg[1]<-paste("<",cc[1])
   my.leg[ncolor]<-paste(">",cc[ncolor-1])
   for (i in 1:(ncolor-2))
     {
      my.leg[i+1]<-paste("[",cc[i],", ",cc[i+1],"[",sep="") 
     }

   legend("bottomright", legend=my.leg, fill=brewer.pal(ncolor,"YlOrBr"), cex=.6, box.lty=0, title=expression(paste("logit(",pi[ij],")",sep="")))


}

