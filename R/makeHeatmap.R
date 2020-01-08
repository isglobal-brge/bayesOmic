#' Creates a Heatmap using coefficients of Bayesian model
#' 
#' @aliases makeHeatmap
#' @param x object of class 'bayesSNPassoc'
#' @param quantiles 
#' @export


makeHeatmap <- function(x, quantiles=c(.2,.5,.6,.8), ...)
 {

   if (!inherits(x, "bayesOmic")) 
     stop("object must be of class 'bayesOmic' - e.g obtained with 'bayesSNPassoc()' or 'bayesOmicAssoc()' functions.")

   
   temp <- x$res.summary$predicted
   predicted <- matrix(temp, ncol=x$N.groups, nrow=x$N.features)
   colnames(predicted) <- x$names.groups
   rownames(predicted) <- x$names.features
   
   ncolor<-length(quantiles)+1
   heatmap(predicted, 
           breaks=c(-1000, quantile(predicted, quantiles), 1000), 
           col=RColorBrewer::brewer.pal(ncolor, "YlOrBr"), cexRow=.6, cexCol=.8, scale="none", 
           margins=c(12,10), ...)
   
   cc <- round(quantile(predicted, quantiles), 4)
   
   my.leg<-rep(NA, ncolor)
   my.leg[1]<-paste("<",cc[1])
   my.leg[ncolor]<-paste(">",cc[ncolor-1])
   for (i in 1:(ncolor-2))
   {
      my.leg[i+1]<-paste("[",cc[i],", ",cc[i+1],"[",sep="") 
   }
   
   legend("bottomright", legend=my.leg, fill=RColorBrewer::brewer.pal(ncolor,"YlOrBr"), 
          cex=.6, box.lty=0, title="Predicted values")
   
   
   
}

