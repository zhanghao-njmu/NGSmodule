args <- commandArgs(TRUE)
rfile <- as.character(args[1])
gfile <- as.character(args[2])
ploid <- as.numeric(args[4])
sample <- as.character(args[5])

library(rtracklayer)
readcount<- import.wig(rfile)
gc<- import.wig(gfile)



rawratios <- read.table(input, header = TRUE)
ratios <- rawratios[order(rawratios$chrmarker), ]
library(DNAcopy)
cols <- names(ratios)
for (i in 3)
{
  # eval
  ratiostmp <- CNA(, ratios$chr, ratios$pos, data.type = "logratio", sampleid = cols[i], presorted = TRUE)
  smooth.ratiostmp <- smooth.CNA(ratiostmp)
  # 	seg.smooth.ratiostmp<-segment(smooth.ratiostmp, verbose=1);
  seg.smooth.ratiostmp <- segment(smooth.ratiostmp, verbose = 1, alpha = 0.001, undo.splits = "sdundo", undo.SD = 1)
  # change alpha for different sensitivity

  png(file = paste(cols[i], "png", sep = "."), width = 1500, height = 480, bg = "white")
  plot(seg.smooth.ratiostmp, ylim = c(-3, 3))
  dev.off()
  write.table(segments.summary(seg.smooth.ratiostmp), quote = FALSE, sep = "\t", file = paste(cols[i], "seg", sep = "."))
}



