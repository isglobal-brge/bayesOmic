print.parameters <-
function(x, ...)
 {
  cat("\n Intercepts (alpha): \n")
  print(x$alpha)
  cat("\n Shared component (log-lambda): \n")
  print(x$loglambda)
 }

