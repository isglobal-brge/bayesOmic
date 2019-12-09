plot.bayesCNVassoc <- function (x, mfrow, pos.legend = "bottomright", ...)
{
    if (!inherits(x, "bayesCNVassoc"))
        stop("object must be of class 'bayesCNVassoc'")

    opar <- par("mfrow", "mar", "xpd")
    on.exit(par(opar))

    N.cnvs <- x$N.cnvs
    N.groups <- x$N.groups
    names.groups <- x$names.groups
    v.median.ind <- x$res.summary
    if (N.groups == 2)
        par(mfrow = c(2, 1))
    else if (N.groups == 3)
        par(mfrow = c(3, 1))
    else if (N.groups == 4)
        par(mfrow = c(2, 2))
    else if (N.groups%in%c(5,6))
        par(mfrow = c(3, 2))
    else {
        stop("more than 6 groups needs 'mfrow' argument to be specified")
        par(mfrow = mfrow)
    }

    par(mar = c(1.5, 2.5, 1, 2))

    ini <- seq(1, nrow(v.median.ind), by=N.cnvs)
    end <- seq(N.cnvs, nrow(v.median.ind), by=N.cnvs)

    conf <- paste("IC", round((1-x$alpha.corrected)*100, 4), "%", sep="")

    mycol <- c("gray80", "green", "red")
    for (i in 1:N.groups) {
        data.i <- v.median.ind[ini[i]:end[i],]
        plot(1:N.cnvs, seq(min(data.i), max(data.i), length.out = N.cnvs), type = "n", xlab = "CNVs", ylab = "", main = names.groups[[i]], axes=FALSE, cex.main=1.4)
        axis(2, cex.axis=1.4)
        col.i <- ifelse(data.i[, 1] > 0, mycol[2], ifelse(data.i[, 3] < 0, mycol[3], mycol[1]))

        for (i in 1:3)
         {
          sel <- col.i == mycol[i]
          points(c(1:N.cnvs)[sel], data.i[sel, 2], pch = 19, cex = 0.9, col=col.i[sel])
          segments(c(1:N.cnvs)[sel], data.i[sel, 1], c(1:N.cnvs)[sel], data.i[sel, 3], col = col.i[sel])
         } 

        segments(0, 0, N.cnvs, 0)
    }
      legend(pos.legend, legend = c(paste(conf, "> 0"), paste(conf, "< 0")),
            lty = c(1, 1), col = mycol[2:3])
      par(xpd=TRUE) 
      pos <- par("usr")
      text(pos[2]/2, pos[3], "CNVs", cex=2.2)
}

