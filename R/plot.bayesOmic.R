#' plot.bayesOmicAssoc
#' 
#' @param x object of class 'bayesOmic'
#' @param mfrow 
#' @param pos.legend position of the plot legend
#' 
#' @S3method plot bayesOmic

plot.bayesOmic <- function (x, type="specific", ...)
{
  if (!inherits(x, "bayesOmic"))
    stop("object must be of class 'bayesOmic'")
  
  type.sel<- charmatch(type, c("specific", "shared"))
  if (is.na(type.sel))
    stop("'type' argument must be 'specific' or 'shared'")
  
  if (attr(x, "control.group"))
    N.groups <- x$N.groups -1
  else
    N.groups <- x$N.groups
  
  if (type.sel==1) {
    lambda <- x$res.summary$lambda
    df <- plyr::ldply (lambda, data.frame) %>% 
      add_column(feature=rep(x$names.features, N.groups))
    names(df)[1:5] <- c("group", "inf", "median", "sup", "sig")
    
    df$mycol <- ifelse(df$sig==0, "gray", ifelse(df$sig=="-1", "red", "blue"))
    plt <- ggplot(df, aes(y=feature, x=median, col=mycol)) + 
      geom_errorbarh(aes(xmin=inf, xmax=sup)) + ylab("Feature") + 
      xlab("Median and Credible Interval") + facet_grid(group ~ .) 
  }
  
  if (type.sel==2){
    df <- data.frame(x$res.summary$u.stats)
    names(df) <- c("mean", "inf", "median", "sup")
    df$feature <- rownames(df)
    df$mycol <- ifelse(df$inf>0, "red", ifelse(df$sup<0, "blue", "lightgray"))
    plt <- ggplot(df, aes(y=feature, x=median, col=mycol)) +
      geom_errorbarh(aes(xmin=inf, xmax=sup)) + ylab("Feature") + 
      xlab("Median and Credible Interval")
      
  }
  
  plt

  }

