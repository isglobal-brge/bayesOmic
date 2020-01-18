library(bayesOmic)
data(armengol)
mod <- bayesOmicAssoc(group="pop", data=armengol)
plot(mod2, type="specific")
plot(mod, type="shared")
x<-mod

library(airway)
data(airway)
airway

roi<-c("2", "100000-1100000")

range <- GRanges(seqnames = roi[1], ranges=roi[2])
airway2<-subsetByOverlaps(airway, range)

mod <- bayesOmicAssoc(group="dex", data=airway2)

mod2 <- bayesOmicAssoc(group="dex", data=airway, roi=roi)

GRanges(airway)
type(airway)
class(airway)

seqnames(rowRanges(airway))
