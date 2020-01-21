library(bayesOmic)
library(GenomicRanges)
library(IRanges)
library(airway)
data(airway)

roi <- GRanges(seqnames = "1", ranges= IRanges(1, 100000))
data<-subsetByOverlaps(airway, roi)
data
mod <- bayesOmicAssoc(group="dex", data=airway, roi=roi)

plot(mod, type="specific")
plot(mod, type="specific")

rleList<-runValue(seqnames(rowRanges(airway)))

chr1 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "1")
chr2 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "2")
chr3 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "3")
chr4 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "4")
chr5 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "5")
chr6 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "6")
chr7 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "7")
chr8 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "8")
chr9 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "9")
chr10 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "10")
chr11 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "11")
chr12 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "12")
chr13 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "13")
chr14 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "14")
chr15 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "15")
chr16 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "16")
chr17 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "17")
chr18 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "18")
chr19 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "19")
chr20 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "20")
chr21 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "21")
chr22 <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "22")
chrX <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "X")
chrY <- subset(airway, subset = runValue(seqnames(rowRanges(airway))) == "Y")

chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY

totals = c(5363, 4047, 3101, 2563, 2859, 2905, 2876, 2386, 2323, 2260, 3208, 2818, 1217, 2244, 2080, 2343, 2903, 1127, 2910, 1317, 736, 1263, 2392, 495)

sum(totals)

# bayesOmicAssocByChr USAGE

models<-bayesOmicAssocByChr(group = "dex", data = airway, chrNames = c("21", "Y"))
models

models$mod_chr_Y
models$mod_chr_21

