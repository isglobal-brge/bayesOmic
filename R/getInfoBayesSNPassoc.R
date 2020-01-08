getInfoBayesSNPassoc <- function(x, Ngroups, Nvar, names.groups, names.SNPs, quantiles)
 {

  nn<-colnames(x[[1]])
  
# intercepts
  o<-grep("alpha",nn)
  alpha<-x[,o]
  tt<-summary(alpha)
  alpha.stats<-cbind(tt[[1]][,1],tt[[2]][,c(1,3,5)])
  rownames(alpha.stats)<-names.groups

# coefficient of shared component
  o<-grep("beta",nn)
  beta<-x[,o]
  tt<-summary(beta)
  beta.stats<-cbind(tt[[1]][,1],tt[[2]][,c(1,3,5)])
  rownames(beta.stats)<-names.groups[-1]


# shared component
  index<-Nvar*(Ngroups-1)
  aux.ini<-seq(1,index,Nvar)
  aux.end<-seq(0,index,Nvar)[-1]

  o<-grep("^u",nn)
  u<-x[,o]
  tt<-summary(u, quantiles=quantiles)
  u.stats<-cbind(tt[[1]][,1],tt[[2]][,c(1,3,5)])
  rownames(u.stats)<-names.SNPs


# specific component
  index<-Nvar*(Ngroups-1)
  aux.ini<-seq(1,index,Nvar)
  aux.end<-seq(0,index,Nvar)[-1]

  o<-grep("^v",nn)
  v<-x[,o]
  tt<-summary(v, quantiles=quantiles)
  v.mean<-matrix(tt[[1]][,1],Nvar,(Ngroups-1))
  v.median<-v.median.ind<-list()
  for (i in 1:(Ngroups-1))
   {
     v.median[[i]]<-tt[[2]][aux.ini[i]:aux.end[i],c(1,3,5)]
     v.median.ind[[i]]<-cbind(v.median[[i]],sig=ifelse(v.median[[i]][,1]>0,1,ifelse(v.median[[i]][,3]<0,-1,0)))
     rownames(v.median.ind[[i]]) <- names.SNPs
   }
  names(v.median.ind) <- names.groups[-1]

# predicted (pi)
  o<-grep("^pi",nn)
  tt<-x[,o]
  pi.mean <- apply(tt[[1]],2,mean)
  
  

  ans<-list(alpha.stats=alpha.stats, 
            beta.stats=beta.stats, 
            lambda=v.median.ind,
            u.stats=u.stats, predicted=pi.mean)


  ans


}

