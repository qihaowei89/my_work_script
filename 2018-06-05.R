#!/usr/bin/env Rscript
# by: qihao
# date : 2018-06-06
# Descriptions:  


options(warn = -1)
package_list <- c("optparse","readxl","magrittr","parallel")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# clean enviroment object
rm(list=ls()) 
options(warn = -1)
#  Load essential packages
suppressMessages({
  library(optparse)
  library(exomeCopy)  #need installed from 'source("https://bioconductor.org/biocLite.R") biocLite("exomeCopy")'
  library(parallel)
  library(magrittr)
})


# 解析命令行
if (TRUE){
  option_list <- list(make_option(c("-B", "--bam"),type="character",help="input bam files"), 
                      make_option(c("-b", "--bed"),type="character",help="target ranges"),
                      make_option(c("-s", "--sample"),type="character",help="input sample"),
                      make_option(c("-o", "--outdir"),type="character",help="output directory"),
                      make_option(c("-t", "--threshold"),type="numeric",default=5,help="min number of each chromosome ranges [default %default]"))
  opts <- parse_args(OptionParser(option_list=option_list))
}

# opts = list()
# opts$bam = "ct_1013_bai_2.target.realigned.rmdup.bam"
# opts$bed = "first.panel.bed"
# opts$sample = "ct_1013_HD786_2"
# opts$outdir = "test"
# opts$threshold = 5

###### Check options 
if(length(opts) < 6){
  cat ("Usage: Rscript CNV.R -i <Input bamFile Directory> -b <Input bedFile> -s <Input samples list> -o <Output Directory> -t <threshold>\n")
  cat ("Example: Rscript CNV.R -i /home/inputdir -b xxx.bed -s list -o /home/outdir -t\n")
  cat ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat ("Version: v1.0\n")
  cat ("Date: 2018-06-06\n")
}


######  control samples bam file 
control_bam_dir <- "."
reference.file <- "./database/human_g1k_v37.fasta"

cat("Step1: Reading target ranges from first.panle.bed\n")
####### read target ranges from first.panle.bed  
target.file <- opts$bed
# target.file <- "first.panel.bed"
target.df   <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end"))
target      <- GRanges(seqname=target.df$seqname,IRanges(start=target.df$start+1,end=target.df$end))
target.sub  <- reduce(target,min.gapwidth=0)

###### Filter chromosome range number less than 5   
target.sub  <- target.sub[target.sub@seqnames %in% target.sub@seqnames@values[target.sub@seqnames@lengths >= 5]]
seqlevels(target.sub) <- target.sub@seqnames@values %>% as.character()

###### Filter  chromosome  with one kind of range width  
seqs_list <- target.sub@seqnames@values
for (i in seqs_list) {
  tmp =target.sub[seqnames(target.sub) == i]
  if (tmp@seqnames@lengths >= 5){
    if ( width(tmp) %>%  unique() %>% length() == 1) {
      start_tmp = tmp@ranges@start
      width_tmp = tmp@ranges@width
      width_tmp[length(width_tmp)] = width_tmp[length(width_tmp)] + 1
      tmp = GRanges(seqnames = i,IRanges(start = start_tmp, width =  width_tmp), strand = tmp@strand)
      target.sub[seqnames(target.sub) == i] = tmp
    }else{
      target.sub[seqnames(target.sub) == i] = tmp
    }
  }
}

cat("Step2: Read samlpe bam files from xxx.target.realigned.rmdup.bam\n")
###### Read samlpe bam files from xxx.target.realigned.rmdup.bam
# samples.list <- read.table(opts$sample,header=F,stringsAsFactors = F) %>% unlist()
# sample_bam <- data.frame(samples="ct_1013_HD786_1",bam="ct_1013_HD786_1.target.realigned.rmdup.bam",stringsAsFactors=F)
sample_bam <- data.frame(samples=opts$sample,bam=opts$bam,stringsAsFactors=F)

####### Read control bam files from xxx.target.realigned.rmdup.bam
control.list <- c("ct_1013_bai_1","ct_1013_bai_2","ct_1013_bai_3","ct_1013_bai_4")
control.bam <- data.frame(samples=control.list,bam=sapply(control.list,function(n) sprintf("%s/%s.target.realigned.rmdup.bam",control_bam_dir,n)),stringsAsFactors=F)

df <- rbind(sample_bam,control.bam)
counts <- target.sub

cat("Step3: Counting reads which in targrt ranges from bam files\n")
####### Counting reads which in targrt ranges from bam files 

for (i in 1:nrow(df)) {
  cat(sprintf("       Starting : %s\n",df$samples[i]))
  s0 <- as.numeric(proc.time()[3])
  mcols(counts)[[df$samples[i]]] <- countBamInGRanges(df$bam[i],target.sub,min.mapq = 30)
  s1 <- as.numeric(proc.time()[3])
  cat(sprintf("       Run time : %f\n",as.numeric(s1-s0)))
}

cat("Step4: Generating background read depth and Calculating GC-content\n")
####### Generating background read depth
counts$bg <- generateBackground(control.list, counts, median)
counts$log.bg <- log(counts$bg + .1)

####### Calculating GC-content 
counts$GC <- getGCcontent(target.sub, reference.file)
counts$GC.sq <- counts$GC^2
counts$width <- width(counts)
counts <- counts[counts$bg > 0,]
# seq="7"
# lapply(seqs_list,function(seq) exomeCopy(counts[seqnames(counts) == seq],opts$sample,X.names=c("log.bg","GC","GC.sq","width"),S=0:4,d=2) %>% copyCountSegments())

###### Plotting CNVs into pdf, Compiling and filtering results, Writeing CNVs results into files
system(sprintf("mkdir -p %s/%s.CNVs.plots",opts$outdir,opts$sample))
cat("Step5: Plotting CNVs into pdf\n       Writing CNVs results into files\n")
out = NULL
for (seq in seqs_list) {
  fit <- exomeCopy(counts[seqnames(counts) == seq],opts$sample,X.names=c("log.bg","GC","GC.sq","width"),S=0:6,d=2) 
  CNV_segment <- copyCountSegments(fit)
  pdf(sprintf("%s/%s.CNVs.plots/%s:chr%s.pdf",opts$outdir,opts$sample,opts$sample,seq),width = 10,height = 6)
  plot(fit, points = T, cols = NULL, show.legend = TRUE,
       main = sprintf("exomeCopy predicted segments\n %s: chr%s",opts$sample,seq),
       xlab = "genomic position", ylab = "normalized read count", xlim = NULL, ylim = NULL, cex = 1, lwd = 3)
  dev.off()
  if(any((CNV_segment$copy.count !=2 & CNV_segment$nranges >= opts$threshold))){
    CNV_segment = CNV_segment[which(CNV_segment$copy.count !=2 & CNV_segment$nranges >= 5)]
    tmp = DataFrame(seqnames = CNV_segment@seqnames,ranges = CNV_segment@ranges, CNV_segment@elementMetadata[,c(1,3)])
    out = rbind(out,tmp)
  } 
}
write.table(out,sprintf("%s/%s.CNV.xls",opts$outdir,opts$sample,".xls"),quote = F,sep = "\t",row.names = F)




