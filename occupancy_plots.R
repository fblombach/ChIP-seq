require(rtracklayer)
require(ggplot2)
require(ggbio)
require(grid)
require(cowplot)

#define coordinates for the genomic region to be plotted
start<-135200
end<-138199

#plotting function
chipPlot<- function(start, end){

#RNA bigwig files are raw coverage, normalisation needs to be applied. CPM values are calculated for mRNA only.
#This is because rRNA depletion in Sulfolobus works poorly
norm_fac<-read.table("data/norm_factors.txt", header=T)
nfac_r1<- subset(norm_fac, sample == "expon1")
nfac_r2 <-subset(norm_fac, sample == "expon2")
genome<-import.gff3("genome_data/Wurtzel.gff")
TSS<-import.bed("genome_data/primaryTSS.bed")
CRISPR_loci<-import.bed("genome_data/crispr_loci.bed")

# set theme for plotting traces
theme_chip<- function (base_size = 11, base_family = "") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.background = element_blank(), 
          panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), 
          plot.margin = unit(c(1,1,0,1), "mm"), complete = TRUE)
}

coord <- IRangesList("gi|15896971|ref|NC_002754.1|" = IRanges::IRanges(start, end))

#data available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141286
#path to directory where data are stored
mydir<-("../GSE141290/")
RNAP1 <- import.bw(paste0(mydir, "GSM4200369_Rpo47_expon_r1.normRatio.bw"), selection = coord, as="NumericList")
RNAP2 <- import.bw(paste0(mydir, "GSM4200370_Rpo47_expon_r2.normRatio.bw"), selection = coord, as="NumericList")
spt45_1 <- import.bw(paste0(mydir, "GSM4200371_Spt45_expon_r1.normRatio.bw"), selection = coord, as="NumericList")
spt45_2 <- import.bw(paste0(mydir, "GSM4200372_Spt45_expon_r2.normRatio.bw"), selection = coord, as="NumericList")
elf1_1 <- import.bw(paste0(mydir, "GSM4200373_Elf1_expon_r1.normRatio.bw"), selection = coord, as="NumericList")
elf1_2 <- import.bw(paste0(mydir, "GSM4200374_Elf1_expon_r2.normRatio.bw"), selection = coord, as="NumericList")
cpsf1_1<- import.bw(paste0(mydir, "GSM4200375_CPSF1_expon_r1.normRatio.bw"), selection = coord, as="NumericList")
cpsf1_2<- import.bw(paste0(mydir, "GSM4200376_CPSF1_expon_r2.normRatio.bw"), selection = coord, as="NumericList")
TFB1_1 <- import.bw(paste0(mydir, "GSM4200377_TFB1_expon_r1.normRatio.bw"), selection = coord, as="NumericList")
TFB1_2 <- import.bw(paste0(mydir, "GSM4200378_TFB1_expon_r2.normRatio.bw"), selection = coord, as="NumericList")
TFEb_1 <- import.bw(paste0(mydir, "GSM4200381_TFEb_expon_r1.normRatio.bw"), selection = coord, as="NumericList")
TFEb_2 <- import.bw(paste0(mydir, "GSM4200382_TFEb_expon_r2.normRatio.bw"), selection = coord, as="NumericList")
TFEa_1 <- import.bw(paste0(mydir, "GSM4200379_TFEa_expon_r1.normRatio.bw"), selection = coord, as="NumericList")
TFEa_2 <- import.bw(paste0(mydir, "GSM4200380_TFEa_expon_r2.normRatio.bw"), selection = coord, as="NumericList")

fw1 <- import.bw(paste0(mydir, "GSM4203662_expon1.plus.bw"), selection = coord, as="NumericList")
fw2 <- import.bw(paste0(mydir, "GSM4203663_expon2.plus.bw"), selection = coord, as="NumericList")
rv1 <- import.bw(paste0(mydir, "GSM4203662_expon1.minus.bw"), selection = coord, as="NumericList")
rv2 <- import.bw(paste0(mydir, "GSM4203663_expon2.minus.bw"), selection = coord, as="NumericList")



#merge replicates and calculate mean
TFEb<- data.frame(r1=data.frame(TFEb_1)[,3], r2=data.frame(TFEb_2)[,3])
TFEb$mean <- rowMeans(TFEb)

TFEa<- data.frame(r1=data.frame(TFEa_1)[,3], r2=data.frame(TFEa_2)[,3])
TFEa$mean <- rowMeans(TFEa)

TFB1<- data.frame(r1=data.frame(TFB1_1)[,3], r2=data.frame(TFB1_2)[,3])
TFB1$mean <- rowMeans(TFB1)

spt45<- data.frame(r1=data.frame(spt45_1)[,3], r2=data.frame(spt45_2)[,3])
spt45$mean <- rowMeans(spt45)

elf1<- data.frame(r1=data.frame(elf1_1)[,3], r2=data.frame(elf1_2)[,3])
elf1$mean <- rowMeans(elf1)

cpsf1<- data.frame(r1=data.frame(cpsf1_1)[,3], r2=data.frame(cpsf1_2)[,3])
cpsf1$mean <- rowMeans(cpsf1)

RNAP<- data.frame(r1=data.frame(RNAP1)[,3], r2=data.frame(RNAP2)[,3])
RNAP$mean <- rowMeans(RNAP)

fw<- data.frame(r1=data.frame(fw1)[,3], r2=data.frame(fw2)[,3])
fw$r1_norm<- fw[,1]*10^6/nfac_r1[,3]/nfac_r1[,4]
fw$r2_norm<- fw[,2]*10^6/nfac_r2[,3]/nfac_r2[,4]
fw$mean <- rowMeans(fw[,3:4])

rv<- data.frame(r1=data.frame(rv1)[,3], r2=data.frame(rv2)[,3])
rv$r1_norm<- rv[,1]*10^6/nfac_r1[,3]/nfac_r1[,4]
rv$r2_norm<- rv[,2]*10^6/nfac_r2[,3]/nfac_r2[,4]
rv$mean <- rowMeans(rv[,3:4])



# genome plots
pGenome<-ggplot() + geom_hline(yintercept=1, size=1) +
  geom_arrowrect(genome, arrow.head.fix=250, fill="dimgrey", colour="dimgrey") +
  geom_rect(data=TSS, show.legend=F, rect.height=0.6) + xlim(start, end) +
  geom_arrowrect(CRISPR_loci, arrow.head.fix=250, 
                 rect.height=0.25) +
  theme_void() +  scale_x_continuous(expand = c(0, 0))
scale_data<- data.frame(x1=start, x2=(start+1000), y1=1, y2=1)
scalebar<- ggplot() + geom_segment(data=scale_data, aes(x=x1, y=y1, xend=x2, yend=y2)) +
  theme_void() + xlim(start,end) + scale_x_continuous(expand = c(0, 0))


# line plots with averaged occupancy and SE ribbon
lm <- ceiling(max(RNAP$r1, RNAP$r2, TFB1$r1, TFB1$r2, TFEb$r1, TFEb$r2, 
                spt45$r1, spt45$r2, elf1$r1, elf1$r2, cpsf1$r1, cpsf1$r2)/20)*20
maxRNA <- max(fw$mean, rv$mean)
lm_RNA <- ceiling(maxRNA/(2*10^(floor(log(maxRNA,10))-1)))*2*10^(floor(log(maxRNA,10))-1)
pRNAP<-ggplot(RNAP, aes(x=c(1:(end-start+1)),y=mean)) + geom_line(colour="#000000", size=0.5) + 
  geom_ribbon(aes(ymin=pmin(r1,r2),ymax=pmax(r1,r2)), alpha=0.5, fill="#000000") + theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(lm), limits = c(1, lm), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1)))  
pSpt45<-ggplot(spt45, aes(x=c(1:(end-start+1)),y=mean)) + geom_line(colour="#453781", size=0.5) + 
  geom_ribbon(aes(ymin=pmin(r1,r2),ymax=pmax(r1,r2)), alpha=0.5, fill="#453781") + theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(lm), limits = c(1, lm), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1)))  
pElf1<-ggplot(elf1, aes(x=c(1:(end-start+1)),y=mean)) + geom_line(colour="#C05F98", size=0.5) + 
  geom_ribbon(aes(ymin=pmin(r1,r2),ymax=pmax(r1,r2)), alpha=0.5, fill="#C05F98") + theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(lm), limits = c(1, lm), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1)))  
pCPSF1<-ggplot(cpsf1, aes(x=c(1:(end-start+1)),y=mean)) + geom_line(colour="#664b4b", size=0.5) + 
  geom_ribbon(aes(ymin=pmin(r1,r2),ymax=pmax(r1,r2)), alpha=0.5, fill="#664b4b") + theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(lm), limits = c(1, lm), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1)))  
pTFEa<-ggplot(TFEa, aes(x=c(1:(end-start+1)),y=mean)) + geom_line(colour="#FFA500", size=0.5) + 
  geom_ribbon(aes(ymin=pmin(r1,r2),ymax=pmax(r1,r2)), alpha=0.5, fill="#FFA500") + theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(lm), limits = c(1, lm), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1)))  
pTFEb<-ggplot(TFEb, aes(x=c(1:(end-start+1)),y=mean)) + geom_line(colour="#FF6347", size=0.5) + 
  geom_ribbon(aes(ymin=pmin(r1,r2),ymax=pmax(r1,r2)), alpha=0.5, fill="#FF6347") + theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(lm), limits = c(1, lm), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1)))  
pTFB1<-ggplot(TFB1, aes(x=c(1:(end-start+1)),y=mean)) + geom_line(colour="#1F968B", size=0.5) + 
  geom_ribbon(aes(ymin=pmin(r1,r2),ymax=pmax(r1,r2)), alpha=0.5, fill="#1F968B") + theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(lm), limits = c(1, lm), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1)))  
pRNAlog<-ggplot(fw, aes(x=c(1:(end-start+1)),y=log(mean,2))) + geom_line(colour="red", size=0.5, alpha=0.5) + 
  geom_ribbon(aes(ymin=log(mean-sd,2),ymax=log(mean+sd,2)), alpha=0.25, fill="red") +
  geom_line(data=rv, aes(x=c(1:(end-start+1)),y=log(mean,2)), colour="blue", size=0.5, alpha=0.5) + 
  geom_ribbon(data=rv, aes(ymin=log(mean-sd,2),ymax=log(mean+sd,2)), alpha=0.25, fill="blue") +
  theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(18), limits = c(0, 18), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1))) 
pRNAfill<- ggplot(fw, aes(x=c(1:(end-start+1)),y=mean)) + geom_line(colour="red", size=0.5) +
  geom_area(alpha=0.5,  fill="red") +
  geom_line(data=rv, aes(x=c(1:(end-start+1)),y=mean), colour="blue", size=0.5, alpha=0.5) + 
  geom_area(data=rv, aes(x=c(1:(end-start+1)), y=mean), alpha=0.5,  fill="blue") +
  theme_chip()  + 
  theme(axis.line.x = element_line(color="black", size = 0.5), 
        axis.line.y = element_line(color="black", size = 0.5)) +
  scale_y_continuous(breaks = c(lm_RNA), limits = c(0, lm_RNA), expand = c(0, 0), position="right") +
  scale_x_continuous(expand = c(0, 0), limits = c(1,(end-start+1))) 

#arrange plots horizontally
plots<-plot_grid(pRNAP, pSpt45, pElf1, pCPSF1, pTFB1, pTFEa, pTFEb, pRNAfill, pGenome, 
                   scalebar, ncol=1, align="v")

}
plots2<-chipPlot(start=start, end=end)
plots2
save_plot(paste0("occupancyPlot.expon.",start,"to", end,".pdf"), 
          plots2,base_height=4, base_width=1+(end-start+1)/3000)
