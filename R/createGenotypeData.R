#' Creates a GenotyeData object with genetic and phenotypic data
#' 
#' @aliases createGenotypeData
#' @param geno object of class 'GdsGenotypeReader' (e.g., PLINK in GDS format)
#' @param pheno a data.frame with phenotypic data
#' @param scanIDcol the column having the identification number in the pheno data.frame (must be the same as in the PLINK data) 
#' @export
#' 
#' 
createGenotypeData <- function(geno, pheno, scanIDcol=1){
  colnames(pheno)[scanIDcol] <- "scanID"
  ids.geno <- getVariable(geno, "sample.id")
  ids.pheno <- pheno$scanID
  if(any(!ids.geno%in%ids.pheno)){
    stop("There are individuals in the geno object that do not have
    phenotypic data. Provide phenotypic data of all individuals in 
    the geno object.")
  }
  
  if (!identical(ids.geno, ids.pheno)){
    rownames(pheno) <- pheno$scanID
    pheno <- pheno[ids.geno, ]
  }
  
  scanAnnot <- ScanAnnotationDataFrame(pheno)
  
  ans <- GenotypeData(geno, scanAnnot = scanAnnot)
  ans
}
