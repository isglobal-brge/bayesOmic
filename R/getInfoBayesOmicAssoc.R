getInfoBayesOmicAssoc<-function(x, N.groups, N.features, names.groups, names.features, quantiles)
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
  
# shared component
  o<-grep("^u",nn)
  u<-x[,o]
  tt<-summary(u, quantiles=quantiles)
  u.stats<-cbind(tt[[1]][,1],tt[[2]][,c(1,3,5)])
  rownames(u.stats)<-names.features  

# predicted
  o<-grep("lambda",nn)
  tt<-x[,o]
  predicted <- apply(tt[[1]],2,mean)


# sigma
  o<-grep("sigma",nn)
  tt<-x[,o]
  sigma<-apply(tt[[1]],2,mean)


# specific component
  index<-N.features*N.groups
  aux.ini<-seq(1,index,N.features)
  aux.end<-seq(0,index,N.features)[-1]

  o<-grep("^v",nn)
  v<-x[,o]
  tt<-summary(v, quantiles=quantiles)
  v.mean<-matrix(tt[[1]][,1],N.features,N.groups)
  v.median<-v.median.ind<-list()
  for (i in 1:N.groups)
   {
     v.median[[i]]<-tt[[2]][aux.ini[i]:aux.end[i],c(1,3,5)]
     v.median.ind[[i]]<-cbind(v.median[[i]],
                              sig=ifelse(v.median[[i]][,1]>0,1,
                                         ifelse(v.median[[i]][,3]<0,-1,0)))
     rownames(v.median.ind[[i]]) <- names.features
  }

  names(v.median.ind) <- names.groups

  ans<-list(alpha.stats=alpha.stats, beta.stats=beta.stats, 
            lambda=v.median.ind, u.stats=u.stats, predicted=predicted, 
            sigma=sigma)

  ans


}


