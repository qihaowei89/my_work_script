#!/usr/bin/env Rscript
# by: qihao
# date : 2018-06-06, 2018-06-27

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


opts = list()
opts$bam = "../control-dp/ct_1013_HD786_1.target.bam"
opts$bed = "/home/wqh/workdir/03_var/ctDNA_103genes_20180702.bed"
opts$sample = "ct_1013_HD786_1"
opts$controlset = "/home/wqh/workdir/03_var/cnv_controlset_ctDNA.RData"
opts$outdir = "./"
opts$threshold = 0


if (TRUE){
  option_list <- list(make_option(c("-B", "--bam"),type="character",help="input bam file"), 
                      make_option(c("-b", "--bed"),type="character",help="target region bed file"),
                      make_option(c("-s", "--sample"),type="character",help="sample name"),
                      make_option(c("-c", "--controlset"),type="character",help="control set"),
                      make_option(c("-o", "--outdir"),type="character",help="output directory"),
                      make_option(c("-t", "--threshold"),type="numeric",default=3,help="copy number threshold to output [default %default]"))
  opts <- parse_args(OptionParser(option_list=option_list))
}

###### Check options 
if(length(opts) < 7){
  cat ("Usage: Rscript cnv_ctDNA.R -B <Input bamFile> -b <Input bedFile> -s <Input samples list> -c <Control set> -o <Output directory> -t <Threshold>\n")
  cat ("Example: Rscript cnv_ctDNA.R -B ct_1013_HD786_1.rmdup.realigned.bam -b first.panel.bed -s ct_1013_HD786_1 -c cnv_controlset_ctDNA.RData -o ./ -t 3\n")
  cat ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat ("Version: v1.0\n")
  cat ("Date: 2018-06-28\n")
}

cat("Step1: Reading target ranges from bed file\n")
####### read target ranges from bed file  
target.file <- opts$bed
target.df   <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end","gene"))
target      <- GRanges(seqname=target.df$seqname,IRanges(start=target.df$start+1,end=target.df$end))
target.sub  <- reduce(target,min.gapwidth=0)
target.gene <- target.sub
mcols(target)$gene <- target.df$gene
target.gene <- target[(target %>% start()) %in% (reduce(target,min.gapwidth=0) %>% start()),]

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

###### Read sample bam file
counts <- target.sub

cat("Step2: Counting reads which in targrt ranges from bam file\n")
####### Counting reads which in targrt ranges from bam file
cat(sprintf("       Starting : %s\n",opts$sample))
s0 <- as.numeric(proc.time()[3])
mcols(counts)[[opts$sample]] <- countBamInGRanges(opts$bam,target.sub,min.mapq = 30)
s1 <- as.numeric(proc.time()[3])
cat(sprintf("       Run time : %f\n",as.numeric(s1-s0)))

cat("Step3: Loading  background read depth and Calculating GC-content from cnv_controlset_ctDNA.RData\n")
control_baseline <- opts$controlset
load(control_baseline)
counts$log.bg <- baseline$log.bg
counts$GC <- baseline$GC 
counts$GC.sq <- baseline$GC.sq
counts$width <- baseline$width
rm(baseline)
#save(counts,file="sample.RData")
getwd()
###### Plotting CNVs into pdf, Compiling and filtering results, Writeing CNVs results into files
system(sprintf("mkdir -p %s/%s.CNVs.plots",opts$outdir,opts$sample))
cat("Step4: Plotting CNVs into pdf\n       Writing CNVs results into files\n")
out = NULL
for (seq in seqs_list) {
  #print(seq)
  range0 = counts[seqnames(counts) == seq]
  range = range0[!(range0@elementMetadata[,1] == 0)]
  if(length(range) >= 5){
    fit <- exomeCopy(range,opts$sample,X.names=c("log.bg","GC","GC.sq","width"),S=0:100,d=2) 
    target0 = target.gene[target.gene@seqnames == seq]
    target = target0[!(range0@elementMetadata[,1] == 0)]
    mcols(target)[["copynum"]] = fit@O.norm
    genes = target@elementMetadata$gene %>% unique() %>% as.character()
    for (gene in genes) {
      tmp = target[target@elementMetadata$gene == gene] %>% sort()
      copy_num = 2*tmp@elementMetadata$copynum %>% median() %>% round(digits = 2) 
      x = apply(ranges(tmp) %>% as.matrix(), 1, function(n) {n[1]+n[2]/2})
      y = 2*tmp@elementMetadata$copynum
      MAX = max(y)
      if(copy_num >= opts$threshold){
        #print(gene)
        out_tmp = DataFrame(chromosome = tmp@seqnames@values,gene = gene,ranges = IRanges(start = start(tmp)[1],end=end(tmp)[length(tmp)]), copy.num = copy_num,nranges = length(tmp))
        out = rbind(out,out_tmp)
        pdf(sprintf("%s/%s.CNVs.plots/%s_%s.pdf",opts$outdir,opts$sample,opts$sample,gene),width = 10,height = 6)
        plot(x,y ,xlab = sprintf("Position of chromosome %s",tmp@seqnames@values),ylab="Copy number", main = sprintf("Copy number of gene  %s",gene),ylim=c(0,MAX+2))
        abline(h = 2, col ="red",lty = 3,type = "a")
        abline(h = copy_num, col ="blue",lty = 4,type = "b")
        dev.off()
      }
    }
  }
}

#reorder output file 
Chromosome=as.character(out$chromosome) %>% as.numeric() 
Chromosome[is.na(Chromosome)] = 23
Out = DataFrame(Chromosome=Chromosome,Start=start(out$ranges),End=end(out$ranges),Gene=out$gene,Segment_length=width(out$ranges),Copy_number=out$copy.num,Bin_count=out$nranges)
Out = Out[order(Out$Chromosome,Out$Start),]
Out$Chromosome[Out$Chromosome == 23] ="X"

write.table(Out,sprintf("%s/%s.CNV.xls",opts$outdir,opts$sample,".xls"),quote = F,sep = "\t",row.names = F)
rm(list = ls())
gc(verbose = F)

