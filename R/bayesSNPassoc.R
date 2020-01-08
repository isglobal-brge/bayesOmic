#' Bayesian model to analyze SNP data
#' 
#' @aliases bayesSNPassoc
#' @aliases plot.bayesOmic
#' @aliases print.bayesOmic
#' @param group vector with one value per individual, stores information about group belonging to (Case/Control, PopA/PopB/PopC...)
#' @param data matrix with the genotype of each SNP (columns) corresponding to each individual (rows)
#' @param sep.allele how alleles are separeted. Argument to be passed through SNPassoc::setupSNP() when SNP data is in a matrix or data.frame. Default is "".
#' @param annotation 
#' @param chr 
#' @param call.rate 
#' @param min.freq 
#' @param sig.level significance level
#' @param n.iter.burn.in 
#' @param n.iter 
#' @param thin 
#' @param n.chain 
#' @param other arguments to be passed trough SNPassoc::setupSNP
#' @export


bayesSNPassoc <- function (group, data, sep.allele="", annotation, chr, call.rate = 0.9, min.freq=0.05, 
                           sig.level = 0.05, n.iter.burn.in = 1000, n.iter = 5000, thin = 10, n.chain = 2,
                            ...)
{

    Y <- group
    if (!is.factor(Y))
        stop("response must be a factor variable")
    if (any(is.na(Y)))
        stop("Dependent variable has missing information")
    
    Nindiv <- length(Y)
    names.groups <- levels(Y)
    Ngroups <- length(names.groups)
    
    
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
    

    nn <- apply(N, 1, sum)
    selec1 <- (nn/Nindiv) >= call.rate
    selec2 <- apply(O/N, 1, max) > min.freq
    N <- N[selec1 & selec2, ]
    O <- O[selec1 & selec2, ]
    Nvar <- nrow(N)
    names.SNPs <- rownames(N)
    N.features <- nrow(N)
    
    
    
    data.JAGS <- list(Ngroups = Ngroups, Nvar = Nvar, O = O, N = N)
    
    initials.JAGS <- list()
    
    for (i in 1:n.chain) {
      alpha <- rnorm(Ngroups, 0, 1)
      logbeta <- rnorm(Ngroups - 1, 0, 1)
      u <- rnorm(Nvar, 0, 1)
      v <- matrix(rnorm(Nvar * (Ngroups - 1), 0, 1), Nvar, Ngroups - 1)
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
    
    alpha.corrected <- (sig.level/2)/N.features
    qq <- c(alpha.corrected, 0.25, 0.5, 0.75, 1-alpha.corrected)
    
    res.summary <- getInfoBayesSNPassoc(res, Ngroups, Nvar, names.groups, names.SNPs, qq)
    
    out <- list(res = res, res.summary = res.summary, model = ans,
            N.groups = Ngroups, N.features = Nvar, names.groups = names.groups,
            names.features = names.SNPs)
    

    attr(out, "control.group") <- TRUE
    class(out) <- "bayesOmic"
    out
}

