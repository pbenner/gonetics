
library(rjson)

## plot densities
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

json1 <- fromJSON(file="ATAC-mm10-forebrain-day11.5.json")$Parameters
json2 <- fromJSON(file="ATAC-mm10-forebrain-day13.5.json")$Parameters
json3 <- fromJSON(file="ATAC-mm10-forebrain-day13.5.normalized.json")$Parameters

pdf("README.pdf", height=4, width=5)
par(mar=c(5,4,1,2))
plot(Y ~ X, json1, type="l", ylab="log density", main="", lty=2)
lines(Y ~ X, json2, type="l", col="darkolivegreen4")
lines(Y ~ X, json3, type="l", col="coral4", lty=2)
legend("topright", legend=c("ATAC Forebrain day 11.5 (reference)", "ATAC Forebrain day 13.5", "ATAC Forebrain day 13.5 (normalized)"), col=c("black", "coral4", "darkolivegreen4"), lty=c(2,1,2), bty="n")
dev.off()

## plot tracks
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------

library(Gviz)
library(GenomicFeatures)
library(org.Mm.eg.db)

## coordinates
## -----------------------------------------------------------------------------

chromosome <- "chr1"

from <- 64680000
to   <- 64750000

## genes
## -----------------------------------------------------------------------------

txdb <- loadDb("/project/wig-data/mouse/genes/TxDb.Mmusculus.UCSC.mm10.refGene.sqlite")

genes <- GeneRegionTrack(txdb, name = "refGene", transcriptAnnotation = "symbol")

symbols <- unlist(mapIds(org.Mm.eg.db, gene(genes), "SYMBOL", "ENTREZID", multiVals = "first"))
symbol(genes) <- symbols[gene(genes)]

## axis
## -----------------------------------------------------------------------------

axis <- GenomeAxisTrack(name=chromosome, background.title = "lightgray", showTitle=TRUE)

## data
## -----------------------------------------------------------------------------

track1  <- DataTrack("ATAC-mm10-forebrain-day11.5.rpm.bw",
                     name = "ATAC Day 11.5", background.title = "brown")
track2  <- DataTrack("ATAC-mm10-forebrain-day13.5.rpm.bw",
                     name = "ATAC Day 13.5", background.title = "brown")
track3  <- DataTrack("ATAC-mm10-forebrain-day13.5.normalized.bw",
                     name = "ATAC Day 13.5 (normalized)", background.title = "brown")


## plot
## -----------------------------------------------------------------------------

pdf("README_browser.pdf", height=6, width=7)
plotTracks(c(axis, track1, track2, track3, genes),
           chromosome = chromosome, from = from, to = to,
           type="l", genome="mm10", ylim=c(0,10))
dev.off()

## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
if (FALSE) {
    ## reduce bin size
    system("../bigWigNil/bigWigNil -v --bin-size=50 /project/wig-data/mouse-embryo/ATAC-mm10-forebrain-day11.5.rpm.bw ATAC-mm10-forebrain-day11.5.rpm.bw")
    system("../bigWigNil/bigWigNil -v --bin-size=50 /project/wig-data/mouse-embryo/ATAC-mm10-forebrain-day13.5.rpm.bw ATAC-mm10-forebrain-day13.5.rpm.bw")

    ## quantile normalize track
    system("./bigWigQuantileNormalize -v ATAC-mm10-forebrain-day11.5.rpm.bw ATAC-mm10-forebrain-day13.5.rpm.bw ATAC-mm10-forebrain-day13.5.normalized.bw")

    ## compute densities
    system("ngstat -v exec ~/Source/ngstat/plugins/nonparametric/nonparametric.so Estimate 100 ATAC-mm10-forebrain-day11.5.rpm.bw ATAC-mm10-forebrain-day11.5.json")
    system("ngstat -v exec ~/Source/ngstat/plugins/nonparametric/nonparametric.so Estimate 100 ATAC-mm10-forebrain-day13.5.rpm.bw ATAC-mm10-forebrain-day13.5.json")
    system("ngstat -v exec ~/Source/ngstat/plugins/nonparametric/nonparametric.so Estimate 100 ATAC-mm10-forebrain-day13.5.normalized.bw ATAC-mm10-forebrain-day13.5.normalized.json")
}
