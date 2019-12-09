getParameters <-
function(x)
 {
  alpha<-x$res.summary$alpha
  loglambda<-x$res.summary$loglambda

  ans<-list(alpha=alpha, loglambda=loglambda)
  
  class(ans)<-"parameters"
  
  ans
 }

