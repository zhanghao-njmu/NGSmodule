#!/usr/bin/env Rscript

###Rscript run-featurecounts.R "/data1/hz/rnacocktail_work/hisat2/SRR970593/alignments.sorted.bam" $gtf $threads $sample 

library(Rsubread)
library(edgeR)
library(Rsamtools)

args <- commandArgs( trailingOnly = TRUE )
nthreads <- args[1]
gtfFile <- args[2]
strandspecific <- args[3]
bamFile <- args[4]
outFilePref <- args[5]

reads_type<-testPairedEndBam(bamFile)

if (reads_type) {
    layout <- "fpkm"
}else{
    layout <- "rpkm"
}

outStatsFilePath  <- paste(outFilePref, '.summary',  sep = '');
outCountsFilePath <- paste(outFilePref, '.count', sep = '');
outFpkmFilePath   <- paste(outFilePref, '.', layout,  sep = '');
outTpmFilePath    <- paste(outFilePref, '.tpm',   sep = '');
outCpmFilePath    <- paste(outFilePref, '.log2CPM',   sep = '');

fCountsList = featureCounts(bamFile, annot.ext=gtfFile, isGTFAnnotationFile=TRUE, nthreads=nthreads, isPairedEnd=reads_type,strandSpecific=strandspecific)
dgeList = DGEList(counts=fCountsList$counts, genes=fCountsList$annotation)
logcpm = cpm(dgeList,log = TRUE,prior.count = 2)
fpkm = rpkm(dgeList, dgeList$genes$Length)
tpm = exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

write(paste("Status",bamFile,sep="\t"),file=outStatsFilePath)
write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE,append=T)

featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts)
colnames(featureCounts) <- c("GeneID",outCountsFilePath)
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

featureCPM = cbind(fCountsList$annotation[,1], logcpm)
colnames(featureCPM) <- c("GeneID",outCpmFilePath)
write.table(featureCPM, outCpmFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

featureFPKM<-cbind(fCountsList$annotation[,1], fpkm)
colnames(featureFPKM) <- c("GeneID",outFpkmFilePath)
write.table(featureFPKM, outFpkmFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

featureTPM<-cbind(fCountsList$annotation[,1], tpm)
colnames(featureTPM) <- c("GeneID",outTpmFilePath)
write.table(featureTPM, outTpmFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
