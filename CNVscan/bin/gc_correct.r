args <- commandArgs(TRUE)

depth.f <- args[1]
sex.f <- args[2]
outfile <- args[3]

rt <- read.table(depth.f,header=TRUE,sep="\t")

# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r120
# GC矫正公式
# Depth_correct = Depth_raw * (m/Mx)
# m是所有exon的中位值
# Mx是与该targe具有相同GC含量的所有target的深度的中位值




shift_sex <- function(df){ # NOT USED
    chr.all <- unique(df$chromosome)
    sex <- NULL
    # 常染色体平均深度
    df.autochr <- df[df$chromosome != "chrX" & df$chromosome != "chrY",]
    mean.cov.autochr <- mean(df.autochr$depth)

    # 检查是否存在chrY
    if ("chrY" %in% chr.all){
        df.chrY <- df[df$chromosome == "chrY",]
        depth.chrY <- df.chrY$depth
        mean.cov.chrY <- mean(depth.chrY)
        if (mean.cov.chrY >= 50){
            sex <- "male"
            print("==========================")
            print(paste("the mean cov of chrY is: ",mean.cov.chrY,sep=""))
            print(paste("the mean cov of autochr is: ",mean.cov.autochr,sep=""))
            print("this sample is assumed as male")
            print("==========================")
        }else{
            sex <- "unknown"
        }
    }

    # 判断是否存在chrX
    if ("chrX" %in% chr.all){
        df.chrX <- df[df$chromosome == "chrX",]
        depth.chrX <- df.chrX$depth
        mean.cov.chrX <- mean(depth.chrX)
        r <- mean.cov.chrX / mean.cov.autochr
        if (r >= 0.4 & r <= 0.6){
            sex <- "male"
            print("==========================")
            print(paste("the mean cov of chrX is: ",mean.cov.chrX,sep=""))
            print(paste("the mean cov of autochr is: ",mean.cov.autochr,sep=""))
            print(paste("mean.cov.chrX / mean.cov.autochr is: ",r,sep=""))
            print(paste("this sample is assumed as male"))
            print("==========================")
        }else{
            sex <- "female"
            print("==========================")
            print(paste("the mean cov of chrX is: ",mean.cov.chrX,sep=""))
            print(paste("the mean cov of autochr is: ",mean.cov.autochr,sep=""))
            print(paste("mean.cov.chrX / mean.cov.autochr is: ",r,sep=""))
            print(paste("this sample is assumed as female"))
            print("==========================")
        }
    }

    if.chrX <- "chrX" %in% chr.all
    if.chrY <- "chrY" %in% chr.all
    if (!if.chrX & !if.chrY){
        sex <- "unknown"
        print("==========================")
        print("this sample can not find chrX or chrY")
        print("this sample sex is unknown")
        print("==========================")
    }
    
    return(sex)
}


# 性别矫正
#this.sex <- shift_sex(rt)

sex.df <- read.table(sex.f,header=FALSE,sep="\t")
this.sex <- sex.df[1,2]
print(paste("infered sex is:",this.sex))

# check chr naming
chrom <- rt[1,1]
if (grepl('chr',chrom)){
    # with chr-prefix
    chr_naming <- 'with_prefix'
}else{
    chr_naming <- 'no_prefix'
}

if (this.sex == "male"){
    if (chr_naming == "with_prefix"){
        df.sexchr <- rt[rt$chromosome == "chrX" | rt$chromosome == "chrY",]
        cov <- df.sexchr$depth # sex chr's cov
        rt[rt$chromosome == "chrX" | rt$chromosome == "chrY",]$depth <- cov * 2
    }else{
        df.sexchr <- rt[rt$chromosome == "X" | rt$chromosome == "Y",]
        cov <- df.sexchr$depth # sex chr's cov
        rt[rt$chromosome == "X" | rt$chromosome == "Y",]$depth <- cov * 2
    }
}


median.depth.all <- median(rt$depth)
print(paste("median depth of all exon is",median.depth.all,sep=": "))

# GC矫正
depth.new <- c()
for (i in 1:nrow(rt)){
    depth.raw <- rt[i,5]
    gc <- rt[i,7]
    same.gc <- rt[rt$gc == gc,]
    med.depth <- median(same.gc$depth)
    if (med.depth == 0){
		dep.new <- depth.raw # 如果med.depth=0,x/0会产生NA/Inf
	}else{
		dep.new <- round(depth.raw * (median.depth.all/med.depth),0) # 四舍五入
	}

    depth.new <- c(depth.new,dep.new)
}

rt$depth_corrected <- depth.new

write.table(rt,outfile,quote=FALSE,sep="\t",row.names=FALSE)






