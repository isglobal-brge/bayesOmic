library(bayesOmic)
data(armengol)
mod <- bayesOmicAssoc(group = "pop", data=armengol)

data(sim.data)
Y <- sim.data$caco
X <- sim.data[, -1]

data("armengol")
Y <- armengol[,1]
X <- armengol[,-1]


mod <- bayesSNPassoc(group = y, data=SNPs, sep="")


library(ggmcmc)

S <-  ggs(mod$res)

Nvar <- 40
Ngroups <- 3
nd <- Nvar * Ngroups

group <- rep(c("A", "B", "C"), each=Nvar)

Alpha <- c(2, 1.9, 2.1)
Beta <- c(0.3, 0.4, 0.2)
Lambda <- rep(0, nd)
Lambda[1:4] <- c(1.2, 1.4, 1.1, -0.8)
Lambda[48:51] <- c(-0.4, 0.4, 0.2, 0.2)
Lambda[95:97] <- c(-0.2, 0.4, 0.9)
O <- cbind(rep(Alpha[1], Nvar) + rep(Beta[1], Nvar) + Lambda[1:40] + rnorm(Nvar, sd=0.01),
           rep(Alpha[2], Nvar) + rep(Beta[2], Nvar) + Lambda[41:80] + rnorm(Nvar, sd=0.01),
           rep(Alpha[3], Nvar) + rep(Beta[3], Nvar) + Lambda[81:120] + rnorm(Nvar, sd=0.01))

library(tidyverse)
library(reshape2)
library(ggplot2)


# plot 1
lambda <- mod$res.summary$lambda

df <- plyr::ldply (lambda, data.frame) %>% 
  add_column(feature=rep(mod$names.features, mod$N.groups))
names(df)[1:5] <- c("group", "inf", "mean", "sup", "sig")
df$mycol <- ifelse(df$sig==0, "lightgray", ifelse(df$sig=="-1", "red", "blue"))


ggplot(df, aes(y=feature, x=mean, col=mycol)) + 
  geom_errorbarh(aes(xmin=inf, xmax=sup)) + facet_grid(group ~ .) 

# plot 2
df <- data.frame(predicted=x$res.summary$predicted,
                 feature=rep(x$names.features, x$N.groups),
                 group=rep(x$names.groups, each=x$N.features)
)

ggplot(df, aes(group, feature, fill=predicted)) + geom_raster() + ylab("Features") +
  xlab("Group") + theme(axis.text.y = element_text(size = 7)) + 
  scale_fill_gradient('predict', limits=c(min(df$predicted), max(df$predicted)),
                      breaks = quantile(df$predicted, breaks),
                      low = "lightblue", high = "darkred") 


# plot 3
temp <- x$res.summary$predicted
predicted <- matrix(temp, ncol=x$N.groups, nrow=x$N.features)
colnames(predicted) <- x$names.groups
rownames(predicted) <- x$names.features

ncolor<-length(breaks)+1
heatmap(predicted, 
        breaks=c(-1000, quantile(predicted, breaks), 1000), 
        col=brewer.pal(ncolor, "YlOrBr"), cexRow=.6, cexCol=.8, scale="none", 
        margins=c(12,10), ...)

cc <- round(quantile(predicted, breaks), 4)

my.leg<-rep(NA, ncolor)
my.leg[1]<-paste("<",cc[1])
my.leg[ncolor]<-paste(">",cc[ncolor-1])
for (i in 1:(ncolor-2))
{
  my.leg[i+1]<-paste("[",cc[i],", ",cc[i+1],"[",sep="") 
}

legend("bottomright", legend=my.leg, fill=RColorBrewer::brewer.pal(ncolor,"YlOrBr"), 
       cex=.6, box.lty=0, title="Predicted values")




a <- 1:Nvar
b <- rep(NA,Nvar)


# general way of writting j1, j2, ..., jNgroups 

matrixJ <- NULL
for (i in 1:Ngroups)
{
  k <- i-1
  jAux <-  c(rep(b, k), a, rep(b, Ngroups - k - 1))
  matrixJ <- cbind(matrixJ, jAux)
}

colnames(matrixJ) <- paste("j", 1:Ngroups, sep="")

alpha <- as.factor(rep(1:Ngroups, each = Nvar))

data <- data.frame(matrixJ, alpha=alpha, y=c(O))

ff.ini <- "y ~ alpha -1  + f(j1, model='iid', constr=TRUE, hyper='logtnormal')"

ff.j <- paste("f(j", 2:Ngroups, " ,copy='j1', fixed=FALSE)", sep="" ,collapse=" + ")

formula.inla <- formula(paste(ff.ini, ff.j, sep=" + "))

res <- inla(formula.inla, family=c("normal"), data = data, 
            control.predictor=list(compute=TRUE), 
            control.inla=list(h=1e-4, strategy = "gaussian", 
                              huge=TRUE, numint.maxfeval= 10000))


lambda.temp <-  data$y-res$summary.linear.predictor[,"mean"] 

lambda <- c(res$summary.random$j1[,"mean"],
            lambda.temp[-c(1:Nvar)])

sd <- sqrt(res$summary.linear.predictor[,"sd"]^2 + 1/res$summary.hyperpar[1,1])
se.lambda <- c(res$summary.random$j1[,"sd"],
               sd[-c(1:Nvar)])

alpha.corrected <- 0.05

inf <- lambda - qnorm(1-(alpha.corrected/2))*se.lambda
sup <- lambda + qnorm(1-(alpha.corrected/2))*se.lambda

lambda.out <- cbind(inf, lambda, sup)
res$summary.fixed
head(lambda.out)
