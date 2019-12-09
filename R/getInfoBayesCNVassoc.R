getInfoBayesCNVassoc<-function(x, N.groups, N.cnvs, names.groups, names.CNVs, quantiles)
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
  rownames(beta.stats)<-names.groups

if (TRUE)
 {
# predicted
  o<-grep("predicted",nn)
  tt<-x[,o]
  predicted<-apply(tt[[1]],2,mean)
 }
 

if (FALSE)
 {
# residuos
  o<-grep("residuos",nn)
  tt<-x[,o]
  residuos<-apply(tt[[1]],2,mean)
}
 


# sigma
  o<-grep("sigma",nn)
  tt<-x[,o]
  sigma<-apply(tt[[1]],2,mean)


# specific component
  index<-N.cnvs*N.groups
  aux.ini<-seq(1,index,N.cnvs)
  aux.end<-seq(0,index,N.cnvs)[-1]

  o<-grep("^v",nn)
  v<-x[,o]
  tt<-summary(v, quantiles=quantiles)
  v.mean<-matrix(tt[[1]][,1],N.cnvs,N.groups)
  v.median<-v.median.ind<-list()
  for (i in 1:N.groups)
   {
     v.median[[i]]<-tt[[2]][aux.ini[i]:aux.end[i],c(1,3,5)]
     v.median.ind[[i]]<-cbind(v.median[[i]],sig=ifelse(v.median[[i]][,1]>0,1,ifelse(v.median[[i]][,3]<0,-1,0)))
   }

# posterior probability of being <> 0
  o<-grep("^p.v",nn)
  temp<-x[,o]
  pv<-pv.stats<-pv.mean<-v.median.select<-v.median.ind.2<-list()
  for (i in 1:N.groups)
   {
    pv[[i]]<-temp[,aux.ini[i]:aux.end[i]]
    pv.mean[[i]]<-apply(pv[[i]][[1]],2,mean)


    v.median.ind.2[[i]]<-cbind(v.median[[i]],"P>1"=pv.mean[[i]])
    colnames(v.median.ind.2[[i]])<-c("2.5%","Median","97.5%","P>0")
    rownames(v.median.ind.2[[i]])<-names.CNVs
   }


  ans<-list(alpha.stats=alpha.stats, beta.stats=beta.stats, v.median.ind=v.median.ind, pv.mean=pv.mean, v.median.ind.2=v.median.ind.2, predicted=predicted, sigma=sigma)

#  ans<-list(alpha.stats=alpha.stats, beta.stats=beta.stats, v.median.ind=v.median.ind, pv.mean=pv.mean, v.median.ind.2=v.median.ind.2, sigma=sigma, residuos=residuos)


  ans


}


