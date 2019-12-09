datChr<-function(x, geno, annotation)
{
  o <- annotation[,1]==x 
  selNames <- annotation[o, 2]
  ans <- geno[, selNames]
  ans
}
