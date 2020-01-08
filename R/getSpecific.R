#' Get specific components 
#' 
#' @param x
#' @param group 
#' 
#' 
#' @export

getSpecific <- function(x, group){
  if (!group%in%x$names.group)
    stop("select a correct group name")
  ss <- x$res.summary$lambda[[group]]
  ans <- ss[ss[,"sig"]!=0,]
  if (nrow(ans)==0){
    cat("No shared components are significant. \n")
    ans <- NULL
  }
  ans
}
  