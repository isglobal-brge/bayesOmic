library(bayesOmic)
library(GenomicRanges)
data(armengol)
mod <- bayesOmicAssoc(group="pop", data=armengol)
plot(mod2, type="specific")
plot(mod, type="shared")
x<-mod

library(airway)
data(airway)
airway

roi1<-c(1, "100000-1100000")
roi2<-c(2, "100000-1100000")
roi<-c("c(1, 2)", "100000-1100000")

range <- GRanges(seqnames = roi[1], ranges=roi[2])
data<-subsetByOverlaps(airway, range)
data

mod1 <- bayesOmicAssoc(group="dex", data=airway, roi=roi1)
mod2 <- bayesOmicAssoc(group="dex", data=airway, roi=roi2)
mod <- bayesOmicAssoc(group="dex", data=airway, roi=roi)

plot(mod1, type="specific")
plot(mod2, type="specific")


GRanges(airway)
type(airway)
class(airway)

seqnames(rowRanges(airway))
