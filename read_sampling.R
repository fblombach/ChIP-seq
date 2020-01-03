# This script samples from all reads to match normal distribution with pre-defined mean and standard deviation
# You will need to convert your input paired-end bam files into bedpe (12 column) format, for example using bamtobedpe with the bedpe option from Bedtools
# The output created is a simple bed file containing the sampled aligments that can be converted back to bam format using bedtobam from Bedtools
# Please note that bamtobed expects bam files to be sorted by name. This can be achieved using "samtools view" using the -n option.
# Likewise you may resort the bam output file using samtools.


# define mean and standard deviation
mean <- 120
sd <- 18

# define name for the output file
name="my_output_file_name.bed"

# read bedpe file generated from bedtools
bedpe <- read.table("my_input_file_name.bedpe")


sampleReads <- function (bedpe, name, mean=120, sd=18){
# Calculate the maximum possible smapling without replacement given the chosen parameters
datalist2=list()
for (i in 51:250){dat2<-nrow(subset(bedpe, V6-V2 == i))  %/%  dnorm(i, mean = mean, sd = sd)
datalist2[[i]] <- dat2
 }
optim = do.call(rbind, datalist2)
size <-  min(optim)
rm(datalist2, dat2)
plot(51:250, optim, log="y", xlab = "fragment size", ylab = "max reads possible to sample", 
     main = name)

#Sampling method: binning data based on fragment size & draw reads from each bin
# for each fragment size i an object dat[i] is added to the datalist
datalist=list()
for (i in 51:250){dat<-bedpe[sample(as.numeric(rownames(subset(bedpe, V6-V2 == i))), size=size*dnorm(i, mean = mean, sd = sd), replace=F),]
datalist[[i]] <- dat
}

# combine the different objects in the data list via rbind to one data frame
bedSamp = do.call(rbind, datalist)
rm(datalist, dat)

# write output as bed file
write.table(bedSamp[,c(1,2,6,7)], file=paste(name, ".mean",mean, ".sd", sd, ".bed", sep=""),quote=F, row.names=F, col.names=F, sep="\t")
print(name)
print(paste("number of read pairs kept"))
print(nrow(bedSamp))
print(paste("number of read pairs in original file"))
print(nrow(bedpe))
rm(bedSamp)
rm(bedpe)
}

sampleReads(bedpe=bedpe, name=name)
