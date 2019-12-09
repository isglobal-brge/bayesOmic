plot.bayesSNPassoc <-
function(x, type="specific", mfrow, ...)
 {

  if (!inherits(x, "bayesSNPassoc")) 
     stop("object must be of class 'bayesSNPassoc'")

  type.sel<- charmatch(type, c("specific", "shared"))
  if (is.na(type.sel))
   stop("'type' argument must be 'specific' or 'shared'")

  Nvar<-x$Nvar

  if (type.sel==1)
   {
    if (attr(x, "control.group"))
     {
      Ngroups<-x$Ngroups - 1
      names.groups<-x$names.groups[-1]
     }
    else
     {
      Ngroups<-x$Ngroups
      names.groups<-x$names.groups
     }

    v.median.ind<-x$res.summary$v.median.ind

    if (Ngroups <= 2)
     par (mfrow=c(2,1)) 
    else if (Ngroups > 2 & Ngroups <= 4)
     par(mfrow=c(2,2))
    else if (Ngroups>4 & Ngroups <= 6)
     par(mfrow=c(3,2))
    else
     {
      stop("more than 6 groups needs 'mfrow' argument to be specified")
      par(mfrow=mfrow)
     }
  
    for (i in 1:Ngroups)
     {
      plot(1:Nvar,seq(min(v.median.ind[[i]][,1:3]),max(v.median.ind[[i]][,1:3]),length.out=Nvar),type="n",xlab="SNPs",ylab="",main=names.groups[[i]])
      points(1:Nvar,v.median.ind[[i]][,1],pch=20,cex=.7)
      points(1:Nvar,v.median.ind[[i]][,2],pch=20,cex=.7)
      points(1:Nvar,v.median.ind[[i]][,3],pch=20,cex=.7)
      segments(1:Nvar,v.median.ind[[i]][,1],1:Nvar,v.median.ind[[i]][,3],col=ifelse(v.median.ind[[i]][,1]>0,3,ifelse(v.median.ind[[i]][,3]<0,4,1))) #Acolorim els significatius
     segments(0,0,Nvar,0)
     legend("bottomright",legend=c("IC95% > 0","IC95% < 0"),lty=c(1,1),col=c(3,4),bty="n",cex=.7)
    }
  }

 if (type.sel==2)
  {
   u.stats<-x$res.summary$u.stats    
   plot(1:Nvar,seq(min(u.stats[,2:4]),max(u.stats[,2:4]), length.out=Nvar), type="n", xlab="SNPs", ylab="", main="Shared component")
   points(1:Nvar,u.stats[,2],pch=20,cex=.7)
   points(1:Nvar,u.stats[,3],pch=20,cex=.7)
   points(1:Nvar,u.stats[,4],pch=20,cex=.7)
   segments(1:Nvar,u.stats[,2], 1:Nvar, u.stats[,4], col=ifelse(u.stats[,2]>0,3,ifelse(u.stats[,4]<0, 4, 1))) #Acolorim els significatius
   segments(0,0,Nvar,0)
   legend("bottomright",legend=c("IC95% > 0","IC95% < 0"),lty=c(1,1),col=c(3,4),bty="n",cex=.7) 
  }

}

