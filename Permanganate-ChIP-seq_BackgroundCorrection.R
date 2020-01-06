# This script performs a background correction for Permanganate ChIP-seq data.
# For each genomic position encoding thymidine (the base with highest sensitivity to permanganate treatment), 
# the n closest neighbouring non-T positions on either side are identied and their median value is subtracted.

# Please note that this script is currently only designed for single chromosomes!
# The expected input format for the raw coverage is bedgraph files split by strand and including coordinates for all nucleotides (lenght of file equals chromosome sequence).

require(seqinr)

data.fw<- "my_permangate_chip-seq_data.fw.bedgraph"
data.rw<- "my_permangate_chip-seq_data.rv.bedgraph"

#Define file with single chromosome sequence
fasta<-read.fasta("my_fasta_file.fasta", forceDNAtolower=T)
genome<-as.vector(unlist(fasta))
genomesize<-length(genome)
genomeComp<- as.data.frame(table(genome))
genomeComp$relFreq<-genomeComp$Freq/sum(genomeComp$Freq)
AT<-mean(genomeComp[c(1,4),3])

## functions used in this script
find_closest_ACG <- function(genome=genome, genomesize=genomesize, n=2) {
  ## return a list with elements corresponding to indices of closest non-T bases
  closest_ACG <- rep(list(NULL), genomesize)
  for (i in 1:genomesize) {
    if (genome[i] == 't') {
      positions <- c()
      j <- i
      while (j > 1 && length(positions) < n) {
        j <- j - 1
        if (genome[j] != 't') positions <- c(positions, j)
      }
      j <- i
      while (j < genomesize && length(positions) < 2*n) {
        j <- j + 1
        if (genome[j] != 't') positions <- c(positions, j)
      }
      closest_ACG[[i]] <- c(positions)
    }
  }
  return(closest_ACG)
}

get_neighbors <- function(data.fw) {
  neighboring_values <- rep(list(NA),genomesize)
  values <- data.fw$V4
  for (i in subset(data.fw, genome == 't')$V3)
    neighboring_values[[i]] <- values[closest_ACG[[i]]]
  return(neighboring_values)
}


# for reverse strand:
find_closest_CGT <- function(genome=genome, genomesize=genomesize, n=2) {
  ## return a list with elements corresponding to indices of closest non-A bases, corresponding to T on reverse strand
  closest_CGT <- rep(list(NULL), genomesize)
  for (i in 1:genomesize) {
    if (genome[i] == 'a') {
      positions <- c()
      j <- i
      while (j > 1 && length(positions) < n) {
        j <- j - 1
        if (genome[j] != 'a') positions <- c(positions, j)
      }
      j <- i
      while (j < genomesize && length(positions) < 2*n) {
        j <- j + 1
        if (genome[j] != 'a') positions <- c(positions, j)
      }
      closest_CGT[[i]] <- c(positions)
    }
  }
  return(closest_CGT)
}

get_neighborsRv <- function(data.rv) {
  neighboring_values <- rep(list(NA),genomesize)
  values <- data.rv$V4
  for (i in subset(data.rv, genome == 'a')$V3)
    neighboring_values[[i]] <- values[closest_CGT[[i]]]
  return(neighboring_values)
}



# find neighbouring non-T positions
closest_ACG <- find_closest_ACG(genome=genome, genomesize=genomesize, n=2) 
closest_CGT <- find_closest_CGT(genome=genome, genomesize=genomesize, n=2) 

# calculate median background
neighboring_values<-get_neighbors(data.fw)
median_nbv<-sapply(neighboring_values, median)

neighboring_valuesRv<-get_neighborsRv(data.rv)
median_nbvRv<-sapply(neighboring_valuesRv, median)

# create bedgraph files
data.fw.corr<-data.frame(data.fw[,1:3], ifelse(genome == "t", data.fw[,4] - median_nbv, 0))
data.rv.corr<-data.frame(data.rv[,1:3], ifelse(genome == "a", data.rv[,4] - median_nbvRv, 0))
# remove negative values (i.e. background is higher than T signal)
data.fw.corr[data.fw.corr[,4] < 0,4]<-0
data.rv.corr[data.rv.corr[,4] < 0,4]<-0

# write out background-corrected data
write.table(data.fw.corr, "my.file.name.fw.bedgraph", 
            sep="\t", quote=F, col.names=F, row.names=F)
write.table(data.rv.corr, "my.file.name.rv.bedgraph", 
            sep="\t", quote=F, col.names=F, row.names=F)

# rescale everything to 1x genome coverage
sfT<-genomesize*AT/(sum(data.fw.corr[,4]) + sum(data.rv.corr[,4]))
data.fw.corr.sc<-data.frame(data.fw[,1:3], data.fw.corr[,4] * sfT)
data.rv.corr.sc<-data.frame(data.rv[,1:3], data.rv.corr[,4] * sfT)

write.table(data.fw.corr.sc, "my.file.name.1xCov.fw.bedgraph", 
            sep="\t", quote=F, col.names=F, row.names=F)

write.table(data.rv.corr.sc, "my.file.name.1xCov.rv.bedgraph", 
            sep="\t", quote=F, col.names=F, row.names=F)

rm(fw, rv, neighboring_values, neighboring_valuesRv, data.fw, data.fw.corr, data.rv, 
   data.rv.corr, median_nbv, median_nbvRv, sfT)  

