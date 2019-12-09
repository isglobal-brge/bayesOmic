snpMatrix2bayes <-
function(x)
 {
   ss<-unclass(x)
   ss<-apply(ss, 2, as.character)
   ss[ss=="00"]<-NA

   ans<-matrix(NA, ncol=ncol(ss), nrow=nrow(ss)) 

   for (i in 1:ncol(ss))
    {
     tt<-table(ss[,i])

     if (length(tt)>1)
      {
       if (tt[1]>tt[length(tt)])
        ans[,i]<-as.numeric(snp(ss[,i], name=c("01","02","03")))-1
       else
        ans[,i]<-as.numeric(snp(ss[,i], name=c("03","02","01")))-1
      }
    }
  colnames(ans)<-colnames(x)
  o<-apply(ans, 2, function(x) sum(is.na(x))==length(x))
  out<-ans[,!o] 
  out
 }

