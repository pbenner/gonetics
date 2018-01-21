
library(rjson)

## plot densities
## -----------------------------------------------------------------------------

json1 <- fromJSON(file="ATAC-mm10-forebrain-day11.5.json")$Distributions[[1]]$Parameters
json2 <- fromJSON(file="ATAC-mm10-forebrain-day13.5.json")$Distributions[[1]]$Parameters
json3 <- fromJSON(file="ATAC-mm10-forebrain-day13.5.normalized.json")$Distributions[[1]]$Parameters

png("README.png", height=350, width=400)
plot(Y ~ X, json1, type="l", ylab="log density", main="")
lines(Y ~ X, json2, type="l", col="coral4")
lines(Y ~ X, json3, type="l", col="darkolivegreen4", lty=2)
legend("topright", legend=c("ATAC Forebrain day 11.5 (reference)", "ATAC Forebrain day 13.5", "ATAC Forebrain day 13.5 (normalized)"), col=c("black", "coral4", "darkolivegreen4"), lty=c(1,1,2), bty="n")
dev.off()

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
