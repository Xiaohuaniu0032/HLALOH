#library(ggplot2)
args <- commandArgs(TRUE)

vaf <- args[1]
name <- args[2]
outdir <- args[3]

rt <- read.table(vaf,header=TRUE,sep="\t")
out.fig <- paste(outdir,'/',name,'.vaf.fig.png',sep="")

ppi <- 80
png(file=out.fig,width=8*ppi,height=44*ppi)
#pdf(file=out.fig,width=6,height=50)
par(mfrow=c(22,1))
par(mar=c(5,5,4,2))
chr.int <- 1:22
for (c in chr.int){
	sub.rt <- rt[rt$chr==c,]
	sub.rt$line <- 1:nrow(sub.rt)
	sub.rt$vaf <- sub.rt$alt_freq
	plot(sub.rt$line,sub.rt$vaf,main=name,xlab=paste("chr",c,sep=""),ylab="vaf",pch=20,cex.lab=2)
}
dev.off()
