#' Bayesian model to analyze CNV data
#' 
#' @aliases bayesOmicAssoc plot.bayesOmic print.bayesOmic
#' @param group name of grouping variable
#' @param data omic data. It can be a data frame, ...
#' @param roi Range of interest. A vector containing the chromosome number and the start-end positions
#' @param miss.indiv minimum percentage of missing individuals allowed. Default is 0.9
#' @param sig.level significance level
#' @param n.iter.burn.in 
#' @param n.iter
#' @param thin
#' @param n.chain
#' @export


bayesOmicAssoc <- function(group, data, roi=NA, miss.indiv = 0.9, sig.level = 0.05, 
                           n.iter.burn.in=1000, n.iter=5000, thin=10,
                           n.chain=2, ...)
 {

   ## Get required data 
   if (is(data, "ExpressionSet")){
      X <- t(Biobase::exprs(data))
      Y <- Biobase::pData(data)[, group]
   } else if (is(data, "SummarizedExperiment")){
     if (!is.na(roi)){
       range <- GRanges(seqnames = roi[1], ranges=roi[2])
       data<-subsetByOverlaps(data, range)
     }
      X <- t(SummarizedExperiment::assay(data))
      Y <- SummarizedExperiment::colData(data)[, group]
   } else if (is.data.frame(data)){
      i <- which(colnames(data)==group)
      X <- data[ , -i]
      Y <- data[ , i]
   } else {
      stop("set must be a data.frame, ExpressionSet or a SummarizedExperiment.")
   }
    
   if (is.null(Y))
      stop("'group' variable must be a name in the data object.")
   if (!is.factor(Y))
      stop("'group' variable must be a factor.")
   
   names.groups<-levels(Y)
   N.groups<-length(names.groups)
   N.features<-ncol(X)
   names.features<-dimnames(X)[[2]]


# Aggregate by group 
      
   data.agg <- t(aggregate(X, by=list(Y), sum, na.rm=TRUE)[,-1])
   N <- t(aggregate(is.na(X) == 0, by = list(Y), sum)[,-1])
   
   O <- data.agg/N
 
# filter features missing individuals
   selec1 <- rep(TRUE, nrow(O))   
   
# filter features by low values (close to 0 ???)
#     selec2<-apply(O/N, 1, max)>0.05
   selec2 <- rep(TRUE, length(selec1))
   O <- O[selec1 & selec2, ] 

   data.JAGS <- list(N.groups=N.groups, N.features=N.features, O=O)
   
   # Initial values
   initials.JAGS<-list()
   for (i in 1:n.chain)
   {
      alpha <- rnorm(N.groups, 0, 1)
      beta <- rnorm(N.groups,0,1)
      u <- rnorm(N.features, 0, 1)
      v <- matrix(rnorm(N.features*N.groups, 0, 1),  N.features, N.groups)
      sigma.v <- abs(rnorm(N.groups,0,1))
      sigma<-rep(1,N.groups)
      initials.JAGS[[i]] <- list(alpha=alpha, beta=beta, sigma=sigma, 
                                 sigma.v=sigma.v, u=u, v=v)
   }
   
   # Run JAGS
   ans<-jags.model(system.file("extdata/JAGSmodels/model_Omic_noControls.bug", package="bayesOmic"), 
                   data=data.JAGS, inits=initials.JAGS, n.chain=n.chain, ...)
   update(ans, n.iter.burn.in)
   
   # Run CODA
   res<-coda.samples(ans, variable.names=c("alpha", "beta", "u", "v", "p.v", 
                                           "lambda","sigma"), 
                     n.iter=n.iter, thin=thin)
   
   alpha.corrected <- (sig.level/2)/N.features
   qq <- c(alpha.corrected, 0.25, 0.5, 0.75, 1-alpha.corrected)
   
   res.summary<-getInfoBayesOmicAssoc(res, N.groups, N.features, 
                                      names.groups, names.features, qq)
   
   out<-list(res=res, N.groups=N.groups, N.features=N.features, names.groups=names.groups,
             names.features=names.features, res.summary=res.summary, 
             data=O, alpha.corrected=alpha.corrected)
   
   attr(out, "control.group") <- FALSE
   class(out)<-"bayesOmic"
   return(out)
 }

