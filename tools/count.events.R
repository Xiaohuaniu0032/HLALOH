library("Rsamtools",quietly = TRUE)
args <- commandArgs(TRUE)

BAMfile <- args[1]
n <- args[2]

x              <- scanBam(BAMfile, index = BAMfile, param=ScanBamParam(what = scanBamWhat(), tag = 'NM'))
#print(x)
readIDs        <- x[[1]][['qname']]
cigar          <- x[[1]][['cigar']]
editDistance   <- unlist(x[[1]][['tag']])
#print(editDistance)
insertionCount <- sapply(cigar, FUN = function(boop) {return(length(grep(pattern = 'I', x = unlist(strsplit(boop, split = '')))))} )
#print(insertionCount)
deletionCount  <- sapply(cigar, FUN = function(boop) {return(length(grep(pattern = 'D', x = unlist(strsplit(boop, split = '')))))} )

indelTotals    <- sapply(cigar, FUN = function(boop) {
	tmp <- unlist(strsplit( gsub("([0-9]+)","~\\1~",boop), "~" ))
	Is  <- grep(pattern = 'I', x = tmp)
	Ds  <- grep(pattern = 'D', x = tmp)
	total <- sum(as.numeric(tmp[(Is-1)])) + sum(as.numeric(tmp[Ds-1]))
	return(total)
})

misMatchCount <- editDistance - indelTotals
#print(misMatchCount)
eventCount <- misMatchCount + insertionCount + deletionCount
names(eventCount) <- 1:length(eventCount)
passed     <- eventCount[which(eventCount <= n)]
y <- readIDs[as.numeric(names(passed))]
#print(y)
#print(table(y))
y <- names(table(y)[which(table(y) == 2)]) # only use paired reads as final useful reads

#print(y) 



