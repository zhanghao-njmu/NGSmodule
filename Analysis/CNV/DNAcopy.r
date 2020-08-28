args = commandArgs(TRUE)
gcfile=as.character(args[1]); 
countfile=as.character(args[2]);
gcfile="/data/database/human/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/2M_CNV/hg38_2M.gc.txt"
countfile="/data/lab/huangmingqian/scCNV/data_P101SC19010617-01-F040-B42-21/1.rawdata/S-Z22-19-9_FKDL190740072-1a-1/S-Z22-19-9_FKDL190740072-1a-1.mapped.rmdup.counts"

gcfile <- read.table(gcfile,sep = "\t",header = F,stringsAsFactors = F)
countfile <- read.table(countfile,sep = "\t",header = F,stringsAsFactors = F)
dat <- merge(x=)


rawratios<-read.table(input, header=TRUE);
ratios<-rawratios[order(rawratios$chrmarker),]
library(DNAcopy);
cols=names(ratios);

for (i in 3 )
{
#eval
	ratiostmp=  CNA( eval (substitute ( cbind( ratios$var), list(var = cols[i] ) ))  , ratios$chr, ratios$pos, data.type="logratio", sampleid=cols[i], presorted=TRUE)  ;
	smooth.ratiostmp<-smooth.CNA(ratiostmp);
#	seg.smooth.ratiostmp<-segment(smooth.ratiostmp, verbose=1);
	seg.smooth.ratiostmp<-segment(smooth.ratiostmp, verbose=1, alpha=0.001,undo.splits="sdundo", undo.SD=1);
# change alpha for different sensitivity

	png(file=paste(cols[i], "png", sep="."), width = 1500, height = 480, bg="white")
	plot(seg.smooth.ratiostmp, ylim=c(-3,3))
	dev.off();
	write.table(segments.summary(seg.smooth.ratiostmp), quote=FALSE, sep="\t", file=paste(cols[i], "seg", sep=".")); 
}
