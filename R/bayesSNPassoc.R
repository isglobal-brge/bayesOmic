bayesSNPassoc <- function (y, snps, annotation, chr, QC = 0.9, min.freq=0.05, method = "JAGS", n.iter.burn.in = 10000, n.iter = 30000, thin = 50, n.chain = 2,  
    ...)
{
    methods <- c("INLA", "JAGS")
    type <- pmatch(toupper(method), methods)
    if (is.na(type))
        stop("Method should be 'INLA' or 'JAGS'")
    Y <- y
    if (!is.factor(Y))
        stop("response must be a factor variable")
    if (any(is.na(Y)))
        stop("Dependent variable has missing information")
    if (inherits(snps, "snp.matrix")) {
        Nindiv <- length(Y)
        names.groups <- levels(Y)
        Ngroups <- length(names.groups)

        if (!missing(chr))
          {
           if (missing(annotation))
            stop("snp.matrix data requires annotation information")
           else
            genosel <- datChr(chr, snps, annotation)
          }
        else
          {
            genosel <- snps
          }

        alleleAgg <- writeFileGenoDat(genosel, allesum = TRUE,
            Y)
        N <- alleleAgg[[3]]
        O <- alleleAgg[[2]]
        nn <- apply(N, 1, sum)
        selec1 <- (nn/Nindiv) >= QC
        selec2 <- apply(O/N, 1, max) > min.freq
        N <- N[selec1 & selec2, ]
        O <- O[selec1 & selec2, ]
        Nvar <- nrow(N)
        names.SNPs <- colnames(genosel)
        names.SNPs <- names.SNPs[selec1 & selec2]
    }
    else {
        o <- is.na(Y)
        if (any(o)) {
            cat("removing individuals with no information in grouping variable ...")
            cat("done \n")
            if (!is.numeric(snps[, 1]))
                X <- data.frame(lapply(snps[!o, ], function(x) as.numeric(x) -
                  1))
            else X <- snps[!o, ]
            Y <- Y[!o]
 }
        else {
            if (!is.numeric(snps[, 1]))
                X <- data.frame(lapply(snps, function(x) as.numeric(x) -
                  1))
            else X <- snps
        }
        names.groups <- levels(Y)
        Ngroups <- length(names.groups)
        names.SNPs <- dimnames(X)[[2]]
        Nindiv <- length(Y)
        Nvar <- ncol(X)
        data.agg <- t(aggregate(X, by = list(Y), sum, na.rm = TRUE))
        data.agg.num <- matrix(as.numeric(data.agg[-1, ]), Nvar,
            Ngroups)
        parameter.n <- t(aggregate(is.na(X) == 0, by = list(Y),
            sum))
        parameter.n <- matrix(as.numeric(parameter.n[-1, ]),
            Nvar, Ngroups)
        colnames(data.agg.num) <- names.groups
        colnames(parameter.n) <- names.groups
        rownames(data.agg.num) <- rownames(data.agg)[-1]
        rownames(parameter.n) <- rownames(data.agg)[-1]
        N <- parameter.n
        O <- data.agg.num
        nn <- apply(N, 1, sum)
        selec1 <- (nn/Nindiv) >= QC
        selec2 <- apply(O/N, 1, max) > min.freq
        N <- N[selec1 & selec2, ]
        O <- O[selec1 & selec2, ]
        Nvar <- nrow(N)
        names.SNPs <- names.SNPs[selec1 & selec2]
    }
    data.JAGS <- list(Ngroups = Ngroups, Nvar = Nvar, O = O,
        N = N)
    if (type == 1) {

        Nvar <- nrow(O) 
        Ngroups <- ncol(O)
        nd <- Nvar * Ngroups
        y.obs <- as.vector(O)
        two.n <- 2 * as.vector(N)
        a <- 1:Nvar
        b <- rep(NA, Nvar)

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

        alpha <- as.factor(rep(1:Ngroups, each = Nvar))
        data <- data.frame(y = y.obs, matrixJ, matrixLambda, alpha = alpha, two.n = two.n)

    
    # general way of writing formula (for any number of groups)

        

     ff.ini <- "y ~ alpha -1 + f(j1, model='iid', constr=TRUE, initial=-1, hyper='logtnormal')"

     ff.j <- paste("f(j", 2:Ngroups, " ,copy='j1', fixed=FALSE)", sep="" ,collapse=" + ")

     ff.lambda <- paste("f(lambda", 1:(Ngroups-1), " ,model='iid', hyper='logtnormal')", sep="" ,collapse=" + ")


     formula.inla <- formula(paste(ff.ini, ff.j, ff.lambda, sep=" + "))

  
     res <- inla(formula.inla, family = "binomial", data = data,
                      Ntrials = two.n, verbose = TRUE,
                      control.inla = list(strategy = "gaussian", huge = FALSE),
                       ...)
 
      out <- list(res = res, model = formula,
            Ngroups = Ngroups, Nvar = Nvar, names.groups = names.groups,
            names.SNPs = names.SNPs)

    }
    if (type == 2) {
        initials.JAGS <- list()
        for (i in 1:n.chain) {
            alpha <- rnorm(Ngroups, 0, 1)
            loglambda <- rnorm(Ngroups - 1, 0, 1)
            u <- rnorm(Nvar, 0, 1)
            v <- matrix(rnorm(Nvar * (Ngroups - 1), 0, 1), Nvar,
                Ngroups - 1)
            sigma.v <- abs(rnorm(Ngroups - 1, 0, 1))
            sigma.lambda <- abs(rnorm(1))
            sigma.u <- abs(rnorm(1))
            initials.JAGS[[i]] <- list(alpha = alpha, loglambda = loglambda,
                u = u, v = v, sigma.u = sigma.u, sigma.v = sigma.v,
                sigma.lambda = sigma.lambda)
        }
        ans <- jags.model(system.file("inst/JAGSmodels/model.bug", package = "bayesOmic"),
            data = data.JAGS, inits = initials.JAGS, n.chain = n.chain,
            ...)
        update(ans, n.iter.burn.in)
        res <- coda.samples(ans, variable.names = c("alpha",
            "loglambda", "u", "v", "pp.v", "pi"), n.iter = n.iter,
            thin = thin)
        res.summary <- getInfoBayesSNPassoc(res, Ngroups, Nvar,
            names.groups, names.SNPs)
        out <- list(res = res, res.summary = res.summary, model = ans,
            Ngroups = Ngroups, Nvar = Nvar, names.groups = names.groups,
            names.SNPs = names.SNPs)
    }
    attr(out, "control.group") <- TRUE
    class(out) <- "bayesSNPassoc"
    out
}

