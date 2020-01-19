library(bayesOmic)
library(GenomicRanges)
data(armengol)
mod <- bayesOmicAssoc(group="pop", data=armengol)
plot(mod1, type="specific")
plot(mod, type="shared")
x<-mod
BiocManager::install("IRanges")
library(airway)
data(airway)
airway

roi <- GRanges(seqnames = "1", ranges= IRranges(1,10000))
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
