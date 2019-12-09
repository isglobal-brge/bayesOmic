getInfoBayesSNPassoc <-
function(x, Ngroups, Nvar, names.groups, names.SNPs)
 {

  nn<-colnames(x[[1]])
  
# intercepts
  o<-grep("alpha",nn)
  alpha<-x[,o]
  tt<-summary(alpha)
  alpha.stats<-cbind(tt[[1]][,1],tt[[2]][,c(1,3,5)])
  rownames(alpha.stats)<-names.groups

# coefficient of shared component
  o<-grep("loglambda",nn)
  loglambda<-x[,o]
  tt<-summary(loglambda)
  loglambda.stats<-cbind(tt[[1]][,1],tt[[2]][,c(1,3,5)])
  rownames(loglambda.stats)<-names.groups[-1]


# shared component
  index<-Nvar*(Ngroups-1)
  aux.ini<-seq(1,index,Nvar)
  aux.end<-seq(0,index,Nvar)[-1]

  o<-grep("^u",nn)
  u<-x[,o]
  tt<-summary(u)
  u.stats<-cbind(tt[[1]][,1],tt[[2]][,c(1,3,5)])
  rownames(u.stats)<-names.SNPs


# specific component
  index<-Nvar*(Ngroups-1)
  aux.ini<-seq(1,index,Nvar)
  aux.end<-seq(0,index,Nvar)[-1]

  o<-grep("^v",nn)
  v<-x[,o]
  tt<-summary(v)
  v.mean<-matrix(tt[[1]][,1],Nvar,(Ngroups-1))
  v.median<-v.median.ind<-list()
  for (i in 1:(Ngroups-1))
   {
     v.median[[i]]<-tt[[2]][aux.ini[i]:aux.end[i],c(1,3,5)]
     v.median.ind[[i]]<-cbind(v.median[[i]],sig=ifelse(v.median[[i]][,1]>0,1,ifelse(v.median[[i]][,3]<0,-1,0)))
   }

# posterior probability of being <> 0
  o<-grep("^pp.v",nn)
  temp<-x[,o]
  pv<-pv.stats<-pv.mean<-v.median.select<-v.median.ind.2<-list()
  for (i in 1:(Ngroups-1))
   {
    pv[[i]]<-temp[,aux.ini[i]:aux.end[i]]
    pv.mean[[i]]<-apply(pv[[i]][[1]],2,mean)


    v.median.ind.2[[i]]<-cbind(v.median[[i]],"P>1"=pv.mean[[i]])
    colnames(v.median.ind.2[[i]])<-c("2.5%","Median","97.5%","P>0")
    rownames(v.median.ind.2[[i]])<-names.SNPs
   }


# pi's
  index<-Nvar*Ngroups
  aux.ini<-seq(1,index,Nvar)
  aux.end<-seq(0,index,Nvar)[-1]

  o<-grep("^pi",nn)
  temp<-x[,o]
  pi<-pi.mean<-list()
  for (i in 1:Ngroups)
   {
    pi[[i]]<-temp[,aux.ini[i]:aux.end[i]]
    pi.mean[[i]]<-apply(pi[[i]][[1]],2,mean)
   }


  ans<-list(alpha.stats=alpha.stats, loglambda.stats=loglambda.stats, u.stats=u.stats, v.median.ind=v.median.ind, pv.mean=pv.mean, v.median.ind.2=v.median.ind.2, pi.mean=pi.mean)


  ans


}

