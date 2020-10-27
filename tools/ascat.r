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

purity_new <- NULL
ploidy_new <- NULL

if (is.null(purity) | is.null(ploidy)) {
	# no solution for purity & ploidy, return NA
	# set purity as 100% and ploidy as 2 when we can not get solution for purity & ploidy
	purity_new <- 1
	ploidy_new <- 2
	print("WARNING:This sample can not get optimal purity & ploidy by ASCAT software")
}else{
	purity_new <- purity
	ploidy_new <- ploidy
}

val <- c(purity_new,ploidy_new)
names(val) <- c("purity","ploidy")

print(paste("purity is:",purity_new,sep=""))
print(paste("ploidy is:",ploidy_new,sep=""))

write.table(val,file=of,quote=FALSE,row.names=TRUE,col.names=FALSE)


