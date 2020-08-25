#!/usr/bin/env Rscript

# setwd("/data/lab/huangmingqian/scCNV/data_P101SC19010617-01-F040-B42-21/")
# vcfFile <-"/data/lab/huangmingqian/scCNV/data_P101SC19010617-01-F040-B42-21/1.rawdata/S-C16d22-1_FKDL190740072-1a-8/GATK3/S-C16d22-1_FKDL190740072-1a-8.bwa.filter.vcf"
# mySample <- "S-22ES-4_FKDL190740072-1a-14"

args = commandArgs(TRUE)
vcfFile=as.character(args[1]);
mySample=as.character(args[2]);

library(vcfR)
library(stringr)
library(dplyr)
library(ggpubr)
vcf <- read.vcfR(vcfFile,verbose = F)

##### Allele balance ##### 
### AD Allelic depths for the ref and alt alleles in the order listed
ad <- extract.gt(vcf, element = 'AD')
gt <- extract.gt(vcf, element = 'GT')
allele1 <- masplit(ad, sort = 1,record = 1)
allele2 <- masplit(ad, sort = 1,record = 2)
ad1 <- allele1 / (allele1 + allele2) 
ad2 <- allele2 / (allele1 + allele2)

pdf (paste(mySample,".AlleleFreq.pdf",sep=""), height = 4, width = 6)
histdata<- hist(ad2[,1], breaks = seq(0,1,by=0.01), col = "#1f78b4", xaxt="n",main=mySample,xlab="Minor/Major allele frequency")
hist(ad1[,1], breaks = seq(0,1,by=0.01), col = "#a6cee3", add = TRUE)
axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))
summary<- table(gt)
summary <- summary[which(str_detect(pattern = "/",string = names(summary)))]
text(x = 0.25,y = max(histdata$counts)*0.7,labels=paste(c("Genetype:",paste0(names(summary),": ",summary)),collapse = "\n"))
dev.off()

gt <- as.data.frame(gt)
gt[,"chr"]<- rownames(gt) %>% gsub(pattern = "_.*",replacement = "",perl = T)
snp_chr<- table(gt) %>% as.data.frame(stringsAsFactors=F) 
het_chr <- snp_chr[which(!snp_chr[[1]] %in% c("1/1","2/2","3/3")),]
het_chr <- aggregate(data=het_chr,Freq~chr,FUN = sum)
het_chr[,"title"] <- "Heterozygous SNP"
het_chr[,"chr"] <- factor(het_chr[,"chr"],levels = paste0("chr",c(1:30,"X","Y")))
p <- ggbarplot(het_chr,x = "chr",y = "Freq",fill = "Freq",label = "Freq")+
  scale_fill_viridis_c(direction = -1)+ylim(0,max(het_chr[,"Freq"]*1.05))+
  labs(x="Chromosome",y="Frequency")+
  facet_grid(~title)+
  theme(panel.border = element_rect(colour = "black",fill = "transparent"),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "right",
        legend.direction = "vertical",
        strip.text.x = element_text(size = 12),
        strip.background = element_rect(fill = "#FCFCFC"))
ggsave(filename = paste(mySample,".HeterFreq.pdf",sep=""),plot = p,height = 6,width = 8)

# ref <- masplit(ad, sort = 0,record = 1)
# alt1 <- masplit(ad, sort = 0,record = 2)
# alt1_freq <- alt1/(alt1+ref)
# ref_freq <- ref/(alt1+ref)
# histdata<- hist(alt1_freq[,1], breaks = seq(0,1,by=0.01), col = "#1f78b4", xaxt="n",main=mySample,xlab="Minor/Major allele frequency")
# hist(ref_freq[,1], breaks = seq(0,1,by=0.01), col = "#a6cee3", add = TRUE)
# 


# ###  focus on the heterozygotes.
# gt <- extract.gt(vcf, element = 'GT')
# hets <- is.het(gt)
# ad[!hets] <- NA
# 
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# ad1 <- allele1 / (allele1 + allele2)
# ad2 <- allele2 / (allele1 + allele2)
# hist(ad2[,1], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n")
# hist(ad1[,1], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
# axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","1/3","3/4",1))
# 
# ### Allele Depth
# ##  filter our data
# ad <- extract.gt(vcf, element = 'AD')
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# #boxplot(allele1, las=3)
# sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.15, 0.95), na.rm=TRUE)
# # Allele 1
# dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[1,])
# #allele1[dp2 < 0] <- NA
# vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
# dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[2,])
# #allele1[dp2 > 0] <- NA
# vcf@gt[,-1][dp2 > 0] <- NA
# # Allele 2
# dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[1,])
# vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
# dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[2,])
# vcf@gt[,-1][dp2 > 0] <- NA
# 
# ### Now we can check out work with another set of box and whisker plots.
# ad <- extract.gt(vcf, element = 'AD')
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# #boxplot(allele1, las=3)
# 
# ### Now we can see if our histogram of allele balance has been cleaned up.
# gt <- extract.gt(vcf, element = 'GT')
# hets <- is_het(gt)
# ad[!hets] <- NA
# 
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# 
# ad1 <- allele1 / (allele1 + allele2)
# ad2 <- allele2 / (allele1 + allele2)
# 
# hist(ad2[,1], breaks = seq(0,1,by=0.02), col = "#1f78b4", xaxt="n", main="zhanhao")
# hist(ad1[,1], breaks = seq(0,1,by=0.02), col = "#a6cee3", add = TRUE)
# axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), labels=c(0,"1/4","1/3","1/2","2/3","3/4",1))

# 
# ##### Find peaks of density #####
# vcf <- read.vcfR(vcfFile,verbose = F)
# # Filter on depth quantiles.
# sums <- apply(allele1, MARGIN=2, quantile, probs=c(0.1, 0.9), na.rm=TRUE)
# # Allele 1
# dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[1,])
# #allele1[dp2 < 0] <- NA
# vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
# dp2 <- sweep(allele1, MARGIN=2, FUN = "-", sums[2,])
# #allele1[dp2 > 0] <- NA
# vcf@gt[,-1][dp2 > 0] <- NA
# # Allele 2
# dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[1,])
# vcf@gt[,-1][ dp2 < 0 & !is.na(vcf@gt[,-1]) ] <- NA
# dp2 <- sweep(allele2, MARGIN=2, FUN = "-", sums[2,])
# vcf@gt[,-1][dp2 > 0] <- NA
# 
# # Censor homozygotes.
# gt <- extract.gt(vcf, element = 'GT')
# hets <- is_het(gt)
# is.na( vcf@gt[,-1][ !hets ] ) <- TRUE
# 
# 
# # Extract allele depths
# ad <- extract.gt(vcf, element = 'AD')
# allele1 <- masplit(ad, record = 1)
# allele2 <- masplit(ad, record = 2)
# ad1 <- allele1 / (allele1 + allele2)
# ad2 <- allele2 / (allele1 + allele2)
# 
# # Parameters
# #winsize <- 1e5
# #
# winsize <- 4e5
# #bin_width <- 0.1
# #bin_width <- 0.05
# #bin_width <- 0.025
# #
# bin_width <- 0.02
# #bin_width <- 0.01
# 
# 
# # Find peaks
# freq1 <- ad1/(ad1+ad2)
# freq2 <- ad2/(ad1+ad2)
# myPeaks1 <- freq_peak(freq1, getPOS(vcf), winsize = winsize, bin_width = bin_width)
# #myCounts1 <- freq_peak(freq1, getPOS(vcf), winsize = winsize, bin_width = bin_width, count = TRUE)
# is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
# myPeaks2 <- freq_peak(freq2, getPOS(vcf), winsize = winsize, bin_width = bin_width, lhs = FALSE)
# #myCounts2 <- freq_peak(freq2, getPOS(vcf), winsize = winsize, bin_width = bin_width, count = TRUE)
# is.na(myPeaks2$peaks[myPeaks2$counts < 20]) <- TRUE
# 
# 
# for(i in 1:4){
#   hist(freq1[ myPeaks1$wins[i,'START_row']:myPeaks1$wins[i,'END_row'], mySample ],
#        breaks = seq(0,1,by=bin_width), xlim=c(0,1), col=8, main = "", xaxt='n')
#   axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1),
#        labels=c(0,'1/4','1/3','1/2','2/3','3/4',1), las=1)
#   abline(v=myPeaks1$peaks[i,mySample], col=2)
#   if(i==2){ title(main=mySample) }
# }
# 
# 
# ##### Visualization #####
# i <- 1
# 
# layout(matrix(1:2, nrow=1), widths = c(4,1))
# par(mar=c(5,4,4,0))
# 
# mySample <- colnames(freq1)[i]
# plot(getPOS(vcf), freq1[,mySample], ylim=c(0,1), type="n", yaxt='n', 
#      main = mySample, xlab = "POS", ylab = "Allele balance")
# axis(side=2, at=c(0,0.25,0.333,0.5,0.666,0.75,1), 
#      labels=c(0,'1/4','1/3','1/2','2/3','3/4',1), las=1)
# abline(h=c(0.25,0.333,0.5,0.666,0.75), col=8)
# points(getPOS(vcf), freq1[,mySample], pch = 20, col= "#A6CEE344")
# points(getPOS(vcf), freq2[,mySample], pch = 20, col= "#1F78B444")
# segments(x0=myPeaks1$wins[,'START_pos'], y0=myPeaks1$peaks[,mySample],
#          x1=myPeaks1$wins[,'END_pos'], lwd=3)
# segments(x0=myPeaks1$wins[,'START_pos'], y0=myPeaks2$peaks[,mySample],
#          x1=myPeaks1$wins[,'END_pos'], lwd=3)
# 
# bp1 <- hist(freq1[,mySample], breaks = seq(0,1,by=bin_width), plot = FALSE)
# bp2 <- hist(freq2[,mySample], breaks = seq(0,1,by=bin_width), plot = FALSE)
# 
# par(mar=c(5,1,4,2))
# barplot(height=bp1$counts, width=0.02,  space = 0, horiz = T, add = FALSE, col="#A6CEE3")
# barplot(height=bp2$counts, width=0.02,  space = 0, horiz = T, add = TRUE, col="#1F78B4")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
