bayesCNVassoc <-
function(y, cnvs, QC=0.9, method="JAGS", 
         n.iter.burn.in=10000, n.iter=30000, thin=50,
         n.chain=2, alpha.corrected, numint.maxfeval=10000, ...)
 {

    methods <- c("INLA", "JAGS")
    type <- pmatch(toupper(method), methods)
    if (is.na(type)) 
        stop("Method should be 'INLA' or 'JAGS'")


# needs some controls
    Y <- y
    if (!is.factor(Y))
     stop("response must be a factor variable")

    if (any(is.na(Y)))
     stop("Dependent variable has missing information")


#
# NOTE: missing information for this case YES!!!!
#
      o <- is.na(Y) 
      if (any(o))
       {
        cat("removing individuals with no information in grouping variable ...")
        cat("done \n")
        X <- cnvs[!o,]
        Y <- Y[!o] 
       }
      else
       {
        X <- cnvs
       }

      names.groups<-levels(Y)
      N.groups<-length(names.groups)
      N.cnvs<-ncol(X)
      names.CNVs<-dimnames(X)[[2]]


# Aggregate by group 
   data.agg <- t(aggregate(X, by=list(Y), sum, na.rm=TRUE))
   data.agg.num <- matrix(as.numeric(data.agg[-1,]), N.cnvs, N.groups)
   parameter.n <- t(aggregate(is.na(X)==0, by=list(Y), sum))
   parameter.n <- matrix(as.numeric(parameter.n[-1,]), N.cnvs, N.groups)
   colnames(data.agg.num) <- names.groups
   colnames(parameter.n) <- names.groups
   rownames(data.agg.num) <- rownames(data.agg)[-1]
   rownames(parameter.n) <- rownames(data.agg)[-1]

# Data for JAGS

     N <- parameter.n
     O <- data.agg.num / N
     Nindiv<-length(Y)


# filter CNVs QC
     nn <- apply(N, 1, sum)
     selec1 <- (nn/Nindiv)>=QC

# filter CNVs MAF (to be supplied)
#     selec2<-apply(O/N,1, max)>0.05
     selec2<-rep(TRUE, length(selec1))
 
     N <- N[selec1 & selec2, ]
     O <- O[selec1 & selec2, ] 
     Nvar <- nrow(N)
     names.CNVs<-names.CNVs[selec1 & selec2] 
    

# Data for JAGS
   predicted <- matrix(NA,N.cnvs,N.groups)
   data.JAGS <- list(N.groups=N.groups, N.cnvs=N.cnvs, O=O, predicted=predicted)


  if (type==1)
   {

#
# Code provided by Havard Rue     (OJO !!!!! 3 grupos !!!!, esto cambia si +!)
#
   
   Nvar <- nrow(O)
   Ngroups <- ncol(O)
   nd <- Nvar * Ngroups

   a <- 1:Nvar
   b <- rep(NA,Nvar)


   # general way of writting j1, j2, ..., jNgroups 
   matrixJ <- NULL
   for (i in 1:Ngroups)
    {
      k <- i-1
      jAux <-  c(rep(b, k), a, rep(b, Ngroups - k - 1))
      matrixJ <- cbind(matrixJ, jAux)
    }
    colnames(matrixJ) <- paste("j", 1:Ngroups, sep="")
        
    alpha <- as.factor(rep(1:Ngroups, each = Nvar))

   
   
   data <- data.frame(matrixJ, alpha=alpha, y=c(O))




   formula <- y ~ alpha - 1 +
          f(j1, model="iid", constr=TRUE,
                hyper="logtnormal") +
          f(j2, copy="j1", fixed=FALSE) +
          f(j3, copy="j1", fixed=FALSE) 
   
         

   res <- inla(formula, family=c("normal"), data = data, 
               control.predictor=list(compute=TRUE), 
               control.inla=list(h=1e-4, strategy = "gaussian", 
                                 huge=TRUE, numint.maxfeval= numint.maxfeval ), ...)


   lambda <-  data$y-res$summary.linear.predictor[,"mean"] 
  
    
   comp1 <- res$summary.linear.predictor[,"sd"]^2
   comp2 <- 1/res$summary.hyperpar[1,1]
   comp3 <- 1/res$summary.hyperpar[1,2]*as.numeric(!is.na(data$j2))
   se.lambda <- sqrt(comp1+comp2+comp3)
   se.lambda <- sqrt(comp1)


   inf <- lambda - qnorm(1-(alpha.corrected/2))*se.lambda
   sup <- lambda + qnorm(1-(alpha.corrected/2))*se.lambda

   res.summary <- cbind(inf, lambda, sup)


   out<-list(res=res, N.groups=N.groups, N.cnvs=N.cnvs, names.groups=names.groups, names.CNVs=names.CNVs, res.summary=res.summary, data=O, alpha.corrected=alpha.corrected)
   class(out)<-"bayesCNVassoc"

   }

 
  if (type==2)
   { 

# Initial values
   initials.JAGS<-list()
   for (i in 1:n.chain)
     {
      alpha <- rnorm(N.groups, 0, 1)
      beta <- rnorm(N.groups,0,1)
      u <- rnorm(N.cnvs, 0, 1)
      v <- matrix(rnorm(N.cnvs*N.groups, 0, 1),  N.cnvs, N.groups)
      sigma.v <- abs(rnorm(N.groups,0,1))
      sigma<-rep(1,N.groups)
      initials.JAGS[[i]] <- list(alpha=alpha, beta=beta, sigma=sigma, 
                                 sigma.v=sigma.v, u=u, v=v)
     }

# Run JAGS
   ans<-jags.model(system.file("inst/JAGSmodels/model2_noN.bug", package="bayesOmic"), 
                   data=data.JAGS, inits=initials.JAGS, n.chain=n.chain, ...)
   update(ans, n.iter.burn.in)

# Run CODA
   res<-coda.samples(ans, variable.names=c("alpha", "beta", "u", "v", "p.v", 
                                           "predicted","sigma"), 
                     n.iter=n.iter, thin=thin)

   if (missing(alpha.corrected))
    qq<-c(0.025, 0.25, 0.5, 0.75, 0.975)
   else
    qq<-c(alpha.corrected[1], 0.25, 0.5, 0.75, alpha.corrected[2])

   res.summary<-getInfoBayesCNVassoc(res, N.groups, N.cnvs, names.groups, names.CNVs, qq)

   aux<-c(O)-res.summary$predicted 
   ss<- rep(res.summary$sigma,each=N.cnvs)
   residuals<-aux/ss 



   out<-list(res=res, N.groups=N.groups, N.cnvs=N.cnvs, names.groups=names.groups,
             names.CNVs=names.CNVs, res.summary=res.summary, 
             residuals=residuals, data=O, alpha.corrected=alpha.corrected)
   
   attr(out, "control.group") <- TRUE
   class(out)<-"bayesCNVassoc"
  } 
   out 

 }

