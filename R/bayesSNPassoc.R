#' Bayesian model to analyze SNP data
#' 
#' @aliases bayesSNPassoc
#' @aliases plot.bayesOmic
#' @aliases print.bayesOmic
#' @param group vector with one value per individual, stores information about group belonging to (Case/Control, PopA/PopB/PopC...).
#' @param data matrix with the genotype of each SNP (columns) corresponding to each individual (rows).
#' @param sep.allele how alleles are separeted. Argument to be passed through SNPassoc::setupSNP() when SNP data is in a matrix or data.frame. Default is "".
#' @param annotation 
#' @param chr 
#' @param call.rate 
#' @param min.freq 
#' @param sig.level significance level.
#' @param method estimating method. It can be INLA (default) or JAGS.
#' @param n.iter.burn.in 
#' @param n.iter 
#' @param thin 
#' @param n.chain 
#' @param ... other arguments to be passed trough INLA::inla
#' @export


bayesSNPassoc <- function (group, data, sep.allele="", annotation, chr, call.rate = 0.9, min.freq=0.05, 
                           sig.level = 0.05, method="inla", n.iter.burn.in = 1000, n.iter = 5000, thin = 10, n.chain = 2,
                            ...)
{

  methods <- c("JAGS", "INLA")
  type <- pmatch(toupper(method), methods)
  if (is.na(type))
    stop("Method should be 'JAGS' or 'INLA'")
  
  Y <- group
  if (!is.factor(Y))
      stop("response must be a factor variable")
  if (any(is.na(Y)))
        stop("Dependent variable has missing information")
    
  #
  # Get aggregated data
  #
    
    if (inherits(data, "snp.matrix")) {
      if (!missing(chr))
        {
        if (missing(annotation))
          stop("snp.matrix data requires annotation information")
        else
          genosel <- datChr(chr, data, annotation)
        }
        else {
          genosel <- data
        }

        alleleAgg <- writeFileGenoDat(genosel, allesum = TRUE, Y)
        N <- alleleAgg[[3]]
        O <- alleleAgg[[2]]
    }
    
    else if (is(data, "matrix") | is(data, "data.frame")){
      genosel <- SNPassoc::setupSNP(data, 1:ncol(data), sep=sep.allele)
      genosel <- mutate_all(genosel, SNPassoc::additive)
     
      O <- t(aggregate(genosel, by=list(Y), sum, na.rm=TRUE)[,-1])
      N <- t(aggregate(is.na(genosel) == 0, by = list(Y), sum)[,-1])
    }
    
    
    #
    # Performs QC
    #
    
    
    miss <- apply(genosel, 2, function(x) mean(!is.na(x)))
    selec1 <- miss >= call.rate
    selec2 <- apply(O/N, 1, max) > min.freq
    N <- N[selec1 & selec2, ]
    O <<- O[selec1 & selec2, ]
    
    
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
      data.JAGS <- list(Ngroups = Ngroups, Nvar = N.features, O = O, N = N)
    
      initials.JAGS <- list()
    
      for (i in 1:n.chain) {
       alpha <- rnorm(Ngroups, 0, 1)
       logbeta <- rnorm(Ngroups - 1, 0, 1)
       u <- rnorm(N.features, 0, 1)
       v <- matrix(rnorm(N.features * (Ngroups - 1), 0, 1), N.features, Ngroups - 1)
       sigma.v <- abs(rnorm(Ngroups - 1, 0, 1))
       sigma.beta <- abs(rnorm(1))
       sigma.u <- abs(rnorm(1))
       initials.JAGS[[i]] <- list(alpha = alpha, logbeta = logbeta,
                u = u, v = v, sigma.u = sigma.u, sigma.v = sigma.v,
                sigma.beta = sigma.beta)
      }
      
      ans <- jags.model(system.file("extdata/JAGSmodels/model_SNPs_controls.bug", package = "bayesOmic"),
            data = data.JAGS, inits = initials.JAGS, n.chain = n.chain,
            ...)
      update(ans, n.iter.burn.in)
      res <- coda.samples(ans, variable.names = c("alpha", "beta", "u", "v", "pp.v", "pi"), 
                        n.iter = n.iter, thin = thin)

      res.summary <- getInfoBayesSNPassoc(res, Ngroups, N.features, names.groups, names.features, qq)
    
      out <- list(res = res, res.summary = res.summary, model = ans,
            N.groups = Ngroups, N.features = N.features, names.groups = names.groups,
            names.features = names.features)
    
    }
    
    
    if (type==2){
      nd <- N.features * Ngroups
      y.obs <- as.vector(O)
      two.n <- 2 * as.vector(N)
      a <- names.features
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
      ff.ini <- "y ~ alpha -1 + j1 + j2 + j3 + j4"
      
      # removed from initial ... It seems it is OK .. also removed from 'formula.inla' 
      # ff.j <- paste("f(j", 2:Ngroups, " ,copy='j1', fixed=FALSE)", sep="" ,collapse=" + ")
      
      ff.lambda <- paste("f(lambda", 1:(Ngroups-1), " , model='iid', hyper='logtnormal')", sep="" ,collapse=" + ")
      
      formula.inla <- formula(paste(ff.ini,  ff.lambda, sep=" + "))
      
      
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
      res.summary$lambda <- lapply(res$summary.random[-1], ff)
      names(res.summary$lambda) <- names.groups[-1]
      res.summary$predicted <- sapply(res$summary.linear.predictor[,5], function(x) exp(x)/(1+exp(x)))
      
      out <- list(res = res, res.summary=res.summary, model = formula.inla,
                  N.groups = Ngroups, N.features = N.features, names.groups = names.groups,
                  names.features = names.features)
      
    }
    
    attr(out, "control.group") <- TRUE
    class(out) <- "bayesOmic"
    out
}

