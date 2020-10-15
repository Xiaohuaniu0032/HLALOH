library(ASCAT)

args <- commandArgs(TRUE)
tvaf <- args[1]
nvaf <- args[2]
tlogR <- args[3]
nlogR <- args[4]
of <- args[5] # outfile

ascat.bc = ascat.loadData(tlogR,tvaf,nlogR,nvaf)
#ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
#ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)

purity <- ascat.output$aberrantcellfraction
names(purity) <- NULL

ploidy <- ascat.output$ploidy
names(ploidy) <- NULL

print(paste("purity is:",purity,sep=""))
print(paste("ploidy is:",ploidy,sep=""))

val <- c(purity,ploidy)
names(val) <- c("purity","ploidy")
write.table(val,file=of,quote=FALSE,row.names=TRUE,col.names=FALSE)


