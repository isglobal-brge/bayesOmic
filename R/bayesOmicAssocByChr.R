#' Function that performs the bayesOmic analysis by chromosome of a RangedSummarizedExperiment
#' 
#' @param group name of grouping variable
#' @param data omic data. It must be a RangedSummarizedExperiment object
#' @param chrList vector list of all seqnames to analyze
#' @export



bayesOmicAssocByChr <- function(group, data, chrList) {
  
  if (!is(data, "RangedSummarizedExperiment")){
    stop("Set must be a RangedSummarizedExperiment.")
  }
  
  modelList <- list()
  
  for (i in chrList){
    modelName = paste("mod_chr",i, sep = "_")
    modelList[[modelName]] <- bayesOmicAssoc(group = group, data = subset(airway, subset = runValue(seqnames(rowRanges(data))) == i))
    
  }
  
  return(modelList)
}
