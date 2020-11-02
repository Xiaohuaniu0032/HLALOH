library(seqinr,quietly=TRUE) # use read.fasta()
library(Biostrings,quietly=TRUE) # use pairwiseAlignment()


args <- commandArgs(TRUE)
workDir <- args[1]

# output files
hla_a_aln <- paste(workDir,'/hla_a_aln.txt',sep="")
hla_b_aln <- paste(workDir,'/hla_b_aln.txt',sep="")
hla_c_aln <- paste(workDir,'/hla_c_aln.txt',sep="")

hlaPath <- paste(workDir,'/hla.result.new',sep="")
#hlaFa <- paste(workDir,'/patient.hla.fa',sep="")
hlaFa <- paste(workDir,'/patient.hlaFasta.fa',sep="")

hlaAlleles <- read.table(hlaPath, sep = '\t', header = FALSE, as.is = TRUE)
hlaAlleles <- unique(sort(hlaAlleles$V1))
print("hla alleles are:")
print(hlaAlleles)

# read hla fasta
hlaFasta <- read.fasta(hlaFa)

############ Func Start ############
PasteVector <- function(v,sep=""){
	vt <- v[1]
	if(length(v) > 1){
		for(g in 2:length(v)){
			vt <- paste(vt,v[g],sep=sep)
		}
	}

	vt <- paste(vt," EnD",sep="")
	out.v <- sub(" EnD","",vt)
	out.v <- sub("NA , ","",out.v)
	out.v <- sub(" , NA","",out.v)
	out.v <- sub(" , NA , "," , ",out.v)
	return(out.v)
}

getMisMatchPositionsPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE){
	seq1aln <- pattern(alignment) # Get the alignment for the first sequence
	seq2aln <- subject(alignment) # Get the alignment for the second sequence
	
	#print("seq1 aln is:")
	#print(seq1aln)

	#print("seq2 aln is:")
	#print(seq2aln)
	
	seq1aln_base = unlist(strsplit(as.character(seq1aln),split=""))
	#print("seq1aln base is:")
	#print(seq1aln_base)
	seq1aln_base_string <- c()
	for (c in seq1aln_base){
		seq1aln_base_string <- paste(seq1aln_base_string,c,sep="")
	}
	print("seq1aln_base_string is:")
	print(seq1aln_base_string)
	
	
	seq2aln_base = unlist(strsplit(as.character(seq2aln),split=""))
	#print("seq2aln base is")
	#print(seq2aln_base)
	seq2aln_base_string <- c()
	for (c in seq2aln_base){
		seq2aln_base_string <- paste(seq2aln_base_string,c,sep="")
	}
	print("seq2aln_base_string is:")
	print(seq2aln_base_string)
	
	# let's check differences across the whole thing
	#gapsSeq1 <- countPattern("-",as.character(seq1aln))
	#seq1alnresidues <- length(unlist(strsplit(as.character(seq1aln),split="")))-gapsSeq1
	
	
	k <- 1
	seq1Positions <- c()
	
	for (char in unlist(strsplit(as.character(seq1aln),split=""))){
		if (char %in% c('C','G','A','T')){
			seq1Positions <- c(seq1Positions,k)
			k <- k + 1
			next
		}

		if (char %in% c('-')){
			seq1Positions <- c(seq1Positions,k)
			next
		}
	}
	
	k <- 1
	seq2Positions <- c()
	for (char in unlist(strsplit(as.character(seq2aln),split=""))){
		if (char %in% c('C','G','A','T')){
			seq2Positions <- c(seq2Positions,k)
			k <- k + 1
			next
		}
		
		if (char %in% c('-')){
			seq2Positions <- c(seq2Positions,k)
			next
		}
	}
	

	diffSeq1 <- seq1Positions[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]
	diffSeq2 <- seq2Positions[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]
	
	diffType1 <- rep(1,length(diffSeq1))
	diffType1[which(unlist(strsplit(as.character(seq1aln),split=""))[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]%in%'-')] <- 2
	
	diffType2 <- rep(1,length(diffSeq2))
	diffType2[which(unlist(strsplit(as.character(seq2aln),split=""))[unlist(strsplit(as.character(seq2aln),split=""))!=unlist(strsplit(as.character(seq1aln),split=""))]%in%'-')] <- 2
	
	# out list
	out <- list()
	out$diffSeq1 <- diffSeq1
	out$diffSeq2 <- diffSeq2
	out$diffType1 <- diffType1
	out$diffType2 <- diffType2
	return(out)
}

############ Func End ############


for (HLA_gene in c('hla_a','hla_b','hla_c')){
	# get a/b/c alleles
	HLA_As <- grep(HLA_gene,hlaAlleles,value=TRUE)
	if(length(HLA_As)<=1){
		next
		# skip gene with two same alleles
		# homo alleles will not get hla_*_aln.txt file (* stands for a/b/c)
	}
	
	HLA_A_type1 <- HLA_As[1]
	HLA_A_type2 <- HLA_As[2]

	HLA_type1Fasta <- hlaFasta[[HLA_A_type1]]
	HLA_type2Fasta <- hlaFasta[[HLA_A_type2]]

	# perform local pairwise alignment
	seq1 <- PasteVector(toupper(HLA_type1Fasta),sep="")
	print(paste("seq1 is:",HLA_A_type1,sep=""))
	print(seq1)	
	seq2 <- PasteVector(toupper(HLA_type2Fasta),sep="")
	print(paste("seq2 is:",HLA_A_type2,sep=""))
	print(seq2)
	
	sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)	
	tmp <- pairwiseAlignment(seq1, seq2, substitutionMatrix = sigma, gapOpening = -2,gapExtension = -4, scoreOnly = FALSE,type='local')
	
	missMatchPositions <- getMisMatchPositionsPairwiseAlignment(tmp,returnlist=TRUE)
	print("misMatchPos info is:")
	print(missMatchPositions)
	
	diffSeq1 <- missMatchPositions$diffSeq1
	diffSeq2 <- missMatchPositions$diffSeq2
	
	diffType1 <- missMatchPositions$diffType1
	diffType2 <- missMatchPositions$diffType2

	alnDF <- data.frame(diffSeq1=diffSeq1,diffSeq2=diffSeq2,diffType1=diffType1,diffType2=diffType2)
	#alnDF$diffSeq1 <- diffSeq1
	#alnDF$diffSeq2 <- diffSeq2
	#alnDF$diffType1 <- diffType1
	#alnDF$diffType2 <- diffType2
	
	colnames(alnDF) <- c(HLA_A_type1,HLA_A_type2,"HLA_A_type1","HLA_A_type2")
	print(alnDF)

	if ("hla_a" %in% HLA_gene){
		# hla a
		write.table(alnDF,file=hla_a_aln,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
	}
	
	if ("hla_b" %in% HLA_gene){
		# hla b
		write.table(alnDF,file=hla_b_aln,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
	}

	if ("hla_c" %in% HLA_gene){
		# hla c
		write.table(alnDF,file=hla_c_aln,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
	}

	if(length(missMatchPositions$diffSeq1) == 0){
		next
	}
	
	if(length(missMatchPositions$diffSeq1) < 5){
		msg <- 'HLA alleles are very similar (fewer than 5 mismatch positions)! Keep that in mind when considering results.'
		warning(msg)
	}
}
	
			

	
