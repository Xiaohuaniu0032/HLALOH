library(ASCAT)

args <- commandArgs(TRUE)
tvaf <- args[1]
nvaf <- args[2]
tlogR <- args[3]
nlogR <- args[4]
of <- args[5] # outfile

ascat.bc = ascat.loadData(tlogR,tvaf,nlogR,nvaf)
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)


