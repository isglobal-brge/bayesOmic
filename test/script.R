library(rjags)
library(bayesOmic)
data(sim.data)


group<-sim.data$caco
SNPs<-sim.data[,-1]
mod <- bayesSNPassoc(group, SNPs, method = "inla")
