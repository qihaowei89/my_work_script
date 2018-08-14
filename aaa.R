library(exomeCopy)
options(width = 70,digits = 4)

#target.file <- "/home/liubei/workdir/handover/ctDNA_pipeline/data_base/first.panel.bed"
#target.df   <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end")) 
#target      <- GRanges(seqname=target.df$seqname,IRanges(start=target.df$start+1,end=target.df$end))
gr <- GRanges(seqname="1",IRanges(start=1,end=249250621))
target  <- subdivideGRanges(gr,subsize = 5000)
target  <- reduce(target, min.gapwidth=0)
target
#target.sub  <- subdivideGRanges(target)


# scanBamHeader(bam.file)[[1]]$targets
# seqlevels(target)
counts <- target
sample.df <- data.frame(samples=paste("sample",1:length(list.files("sample_set",pattern = ".+bam$")),sep = ""),bam.files=list.files("sample_set",pattern = ".+bam$",full.names = T),stringsAsFactors=FALSE)
contrl.df <- data.frame(samples=paste("con",   1:length(list.files("control_set",pattern = ".+bam$")),sep = ""),bam.files=list.files("control_set",pattern = ".+bam$",full.names = T),stringsAsFactors = F)

df <- rbind(sample.df,contrl.df)
for (i in 1:nrow(df)) {
  mcols(counts)[[df$samples[i]]] <- countBamInGRanges(df$bam.files[i],target,min.mapq = 1)
}


reference.file <- "/home/xiaoshan/Desktop/workdir/database/human_g1k_v37.fasta"
counts$GC <- getGCcontent(target, reference.file)
counts

control.samples <- grep("con.+",colnames(mcols(counts)),value=TRUE)
counts$bg <- generateBackground(control.samples, counts, median)
counts$log.bg <- log(counts$bg + .1)
counts$bg.var <- generateBackground(control.samples, counts, var)
counts
summary(counts$bg)
counts <- counts[counts$bg > 0,]

counts$GC.sq <- counts$GC^2
counts$width <- width(counts)

fit <- exomeCopy(counts[seqnames(counts) == "1"],sample.name="con1",X.names=c("log.bg","GC","GC.sq","width"),S=0:6,d=2)
show(fit)

copyCountSegments(fit)
cnv.cols <- c("red","orange","black","deepskyblue","blue","blue2","blue4")
plot(fit,col=cnv.cols)


vignette("exomeCopy")
