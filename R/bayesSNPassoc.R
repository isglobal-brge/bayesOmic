#' Bayesian model to analyze SNP data
#' 
#' @aliases bayesSNPassoc
#' @aliases plot.bayesOmic
#' @aliases print.bayesOmic
#' @param formula
#' @param data setupSNP data.frame (SNPassoc) with genotypes of each SNP (columns) corresponding to each individual (rows) and phenotypic variables.
#' @param chr 
#' @param call.rate 
#' @param min.freq 
#' @param sig.level significance level.
#' @param method estimating method. It can be INLA (default) or JAGS.
#' @param n.iter.burn.in 
#' @param n.adapt
#' @param n.iter 
#' @param thin 
#' @param n.chain 
#' @param n.cores
#' @param ... other arguments to be passed trough INLA::inla
#' @export


bayesSNPassoc <- function (formula, data,  chr, 
                           call.rate = 0.9, min.freq=0.05, 
                           sig.level = 0.05, method="JAGS", 
                           n.iter.burn.in = 1000, n.adapt=100,
                           n.iter = 5000, thin = 5, n.chain = 2,
                           n.cores = 1, ...)
{

  methods <- c("JAGS", "INLA")
  type <- pmatch(toupper(method), methods)
  if (is.na(type))
    stop("Method should be 'JAGS' or 'INLA'")
  
  group <- all.vars(formula)
    
  #
  # Get aggregated data
  #
    
    if (inherits(data, "GenotypeData")) {
      if (!chr%in%paste(1:24))
        stop("Chromosome must be a character value 1:24")
      mask <- getChromosome(data) == chr
      ss <- which(mask)
      genosel <- getGenotype(data, snp=c(ss[1], length(ss)),
                             transpose=TRUE)
      colnames(genosel) <- getSnpID(data)[mask]
      
      Y <- getScanVariable(data, group)
    }
    
    else if (is(data, "matrix") | is(data, "data.frame")){
      genosel <- data %>% select(labels(data)) %>%
                                   mutate_all(SNPassoc::additive)
      if (!group%in%colnames(data))
        stop(paste(group, "is not in the data object"))
      Y <- data[ , which(colnames(data)==group)]
    }
  else{
   stop("data must be a 'GenotypeData' (GWATools) or a matrix/data.frame (SNPassoc)")
  }
  
  if (!is.factor(Y))
    Y <- as.factor(Y)

  nas <- is.na(Y)
  if (any(nas)) {
    warning("Dependent variable has missing information. Complete cases are analyzed")
    genosel <- genosel[!nas, ]
    Y <- Y[!nas]
  }
    
  O <- t(aggregate(genosel, by=list(Y), sum, na.rm=TRUE)[,-1])
  N <- t(aggregate(is.na(genosel) == 0, by = list(Y), sum)[,-1])
    
    #
    # Performs QC
    #
    
    
    miss <- apply(genosel, 2, function(x) mean(!is.na(x)))
    selec1 <- miss >= call.rate
    selec2 <- apply(O/N, 1, max) > min.freq
    N <- N[selec1 & selec2, ]
    O <- O[selec1 & selec2, ]
    
    
    #
    # Get dimensions
    #
    
    Nindiv <- length(Y)
    names.groups <- levels(Y)
    Ngroups <- length(names.groups)
    names.features <- rownames(N)
    N.features <- nrow(N)
    
    #
    # Correct for multiple comparisons
    #
    
    alpha.corrected <- (sig.level/2)/N.features
    qq <- c(alpha.corrected, 0.25, 0.5, 0.75, 1-alpha.corrected)
    
    if (type==1) {
      
      setts <- list('n.iter' = n.iter, 'n.thin' = thin, 'n.burn' = n.iter.burn.in,
                    'n.chains' = n.chain, 'n.adapt' = n.adapt)
      
      data.JAGS <- list(Ngroups = Ngroups, Nvar = N.features, O = O, N = N)
      
      sharedModel <- system.file("extdata/JAGSmodels/model_SNPs_controls.bug", package = "bayesOmic")           
      params <- c("alpha", "beta", "u", "v", "pi")
      
      if (n.cores>1){
        cl <- makePSOCKcluster(n.cores)
        tmp <- clusterEvalQ(cl, library(dclone))
        res <- jags.parfit(cl = cl, data = data.JAGS,
                         params = params, 
                         model = sharedModel, 
                         n.chains = setts$n.chains, 
                         n.adapt = setts$n.adapt, 
                         n.update = setts$n.burn,
                         n.iter = setts$n.iter, 
                         thin = setts$n.thin)
        stopCluster(cl)
      }
      else{
        res <- jags.fit(data = data.JAGS,
                        params = params, 
                        model = sharedModel, 
                        n.chains = setts$n.chains, 
                        n.adapt = setts$n.adapt, 
                        n.update = setts$n.burn,
                        n.iter = setts$n.iter, 
                        thin = setts$n.thin)
      }
      
      res.summary <- getInfoBayesSNPassoc(res, Ngroups, N.features, 
                                          names.groups, names.features, qq)
    
      out <- list(res = res, res.summary = res.summary, 
            N.groups = Ngroups, N.features = N.features, names.groups = names.groups,
            names.features = names.features)
    
    }
    
    
    if (type==2){
      nd <- N.features * Ngroups
      y.obs <- as.vector(O)
      two.n <- 2 * as.vector(N)
      a <- 1:N.features
      #a <- names.features
      b <- rep(NA, N.features)
      
      # general way of writting j1, j2, ..., jNgroups 
      matrixJ <- NULL
      for (i in 1:Ngroups)
      {
        k <- i-1
        jAux <-  c(rep(b, k), a, rep(b, Ngroups - k - 1))
        matrixJ <- cbind(matrixJ, jAux)
      }
      colnames(matrixJ) <- paste("j", 1:Ngroups, sep="")
      
      matrixLambda <- matrixJ[,-1,drop=FALSE]
      
      lambdas <- paste("lambda", 1:(Ngroups-1), sep="")
      colnames(matrixLambda) <- lambdas
      
      alpha <- as.factor(rep(names.groups, each = N.features))
      data <- data.frame(y = y.obs, matrixJ, matrixLambda, alpha = alpha, two.n = two.n)
      
      
      # general way of writing formula (for any number of groups)
      ff.ini <- "y ~ alpha -1 + f(j1, model='iid', constr=TRUE, initial=-1, hyper='logtnormal')"
      #fs <- paste(paste0("j", 2:Ngroups), collapse="+")
      #ff.ini <- paste("y ~ alpha -1 +", fs)
      
      # removed from initial ... It seems it is OK .. also removed from 'formula.inla' 
      ff.j <- paste("f(j", 2:Ngroups, " ,copy='j1', fixed=FALSE)", sep="" ,collapse=" + ")
      
      ff.lambda <- paste("f(lambda", 1:(Ngroups-1), " , model='iid', hyper='logtnormal')", sep="" ,collapse=" + ")
      
      formula.inla <- formula(paste(ff.ini,  ff.j, ff.lambda, sep=" + "))
      
      
      res <- inla(formula.inla, family = "binomial", data = data,
                  Ntrials = two.n, verbose = FALSE, quantiles=qq,
                  control.predictor=list(compute=TRUE),
                  control.inla = list(strategy = "gaussian", huge = FALSE),
                  keep=FALSE, ...)
      
      res.summary <- list()
      idx <- c(1, 3, 5, 7)
      res.summary$alpha <- res$summary.fixed[, idx]
      
      ff <- function(x, idx=idx2, nn=names.features){
        ans <- x[,idx]
        ans$sig <- ifelse(ans[3] < 0, -1, ifelse(ans[1]>0, 1, 0))
        colnames(ans)[4] <- "sig"
        rownames(ans) <- names.features
        return(ans)
      }
      
      
      idx2 <- c(4, 6, 8)
      res.summary$u.stats <- res$summary.random$j1[,idx+1]
      rownames(res.summary$u.stats) <- res$summary.random$j1$ID
      res.summary$u.stats$sig <- ifelse(res.summary$u.stats[4] < 0, -1, ifelse(res.summary$u.stats[2]>0, 1, 0))
      
      
      nn <- names(res$summary.random)
      ii <- grep("lambda", nn)
      res.summary$lambda <- lapply(res$summary.random[ii], ff)
  
      names(res.summary$lambda) <- names.groups[-1]
      res.summary$predicted <- sapply(res$summary.linear.predictor[,5], function(x) exp(x)/(1+exp(x)))
      
      out <- list(res.summary=res.summary, model = formula.inla,
                  N.groups = Ngroups, N.features = N.features, names.groups = names.groups,
                  names.features = names.features)
      
    }
    
    attr(out, "control.group") <- TRUE
    class(out) <- "bayesOmic"
    out
}

