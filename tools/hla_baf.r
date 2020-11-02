args <- commandArgs(TRUE)
hla_a_baf <- args[1]
hla_b_baf <- args[2]
hla_c_baf <- args[3]
name <- args[4]
outdir <- args[5]

outfig <- paste(outdir,'/',name,'.HLA.BAF.pdf',sep="")

pdf(file=outfig,width=8,height=4)

a <- read.table(hla_a_baf,header=TRUE,sep="\t")
b <- read.table(hla_b_baf,header=TRUE,sep="\t")
c <- read.table(hla_c_baf,header=TRUE,sep="\t")


# plot hla a
if (nrow(a) != 0){
    # have enough het snp pos
    plot(NULL,xlim=c(0,max(a$pos1)),ylim=c(0,1),xlab="Pos Index",ylab="BAF",xaxt="n",xaxs="i",main="HLA_A BAF")
    for (i in 1:nrow(a)){
        val <- a[i,]
        points(val$pos1,val$baf1,pch=20)
        points(val$pos1,val$baf2,pch=20)
    }
}else{
    plot(NULL,xlim=c(0,1000),ylim=c(0,1),xlab="Pos Index",ylab="BAF",xaxt="n",xaxs="i",main="HLA_A BAF")
}
abline(h=0.5,col="blue")


# plot hla b
if (nrow(b) != 0){
    plot(NULL,xlim=c(0,max(b$pos1)),ylim=c(0,1),xlab="Pos Index",ylab="BAF",xaxt="n",xaxs="i",main="HLA_B BAF")
    for (i in 1:nrow(b)){
        val <- b[i,]
        points(val$pos1,val$baf1,pch=20)
        points(val$pos1,val$baf2,pch=20)
    }
}else{
    plot(NULL,xlim=c(0,1000),ylim=c(0,1),xlab="Pos Index",ylab="BAF",xaxt="n",xaxs="i",main="HLA_B BAF")
}
abline(h=0.5,col="blue")

# plot hla c
if (nrow(c) != 0){
    plot(NULL,xlim=c(0,max(c$pos1)),ylim=c(0,1),xlab="Pos Index",ylab="BAF",xaxt="n",xaxs="i",main="HLA_C BAF")
    for (i in 1:nrow(c)){
        val <- c[i,]
        points(val$pos1,val$baf1,pch=20)
        points(val$pos1,val$baf2,pch=20)
    }
}else{
    plot(NULL,xlim=c(0,1000),ylim=c(0,1),xlab="Pos Index",ylab="BAF",xaxt="n",xaxs="i",main="HLA_C BAF")
}
abline(h=0.5,col="blue")


dev.off()
