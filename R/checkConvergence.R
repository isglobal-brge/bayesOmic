checkConvergence <-
function(x, type="Markov chain", parameter="alpha", ...)
 {
  if (!class(x)%in%c("bayesSNPassoc", "bayesCNVassoc")) 
     stop("object must be of class 'bayesSNPassoc'")

  type.plot<-charmatch(type, c("Markov chain", "Gelman-Rubin"))

  if (is.na(type.plot))
   stop("'type' argument must be 'Markov chain' or 'Gelman-Rubin'")


  param<-charmatch(parameter, c("alpha", "log-lambda", "shared", "specific"))
 
  if (is.na(param))
   stop("'parameter' argument must be 'alpha', 'log-lambda', 'shared' or 'specific'")

  nn<-colnames(x$res[[1]])

  if (type.plot==1)
   {
    if (param==1)
     {
      o<-grep("alpha", nn)
      plot(x$res[,o])
     }

    if (param==2)
     {
      o<-grep("loglambda", nn)
      plot(x$res[,o])
     }

    if (param==3)
     {
      o<-grep("^u", nn)
      plot(x$res[,o])
     }

    if (param==4)
     {
      o<-grep("^v", nn)
      plot(x$res[,o])
     }
   }


  if (type.plot==2)
   {
    if (param==1)
     {
      o<-grep("alpha", nn)
      gelman.plot(x$res[,o])
     }

    if (param==2)
     {
      o<-grep("loglambda", nn)
      gelman.plot(x$res[,o])
     }

    if (param==3)
     {
      o<-grep("^u", nn)
      gelman.plot(x$res[,o])
     }

    if (param==4)
     {
      o<-grep("^v", nn)
      gelman.plot(x$res[,o])
     }
   }


 }

