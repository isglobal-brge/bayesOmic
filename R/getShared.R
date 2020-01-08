#' Get shared components 
#' 
#' @param x 
#' 
#' 
#' @export

getShared <- function(x){
  ss <- x$res.summary$u
  sig <- apply(ss, 1, function(x) ifelse(x[2]>0, 1, ifelse(x[4]<0, -1, 0)))
  ans <- ss[sig!=0,]
  if (nrow(ans)==0){
    cat("No shared components are significant. \n")
    ans <- NULL
  }
  ans
}
  