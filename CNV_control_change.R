#!/usr/bin/env Rscript
# by: qihao
# date : 2018-06-06

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
if (T){
  option_list <- list(make_option(c("-B", "--bam"),type="character",help="control bam file list"), 
                      make_option(c("-b", "--bed"),type="character",help="target ranges"),
                      # make_option(c("-c", "--control"),type="character",help="controls' ID"),
                      make_option(c("-o", "--outdir"),type="character",help="output directory")
                      )
  opts <- parse_args(OptionParser(option_list=option_list))
}

opts = list()
opts$bam = "bam_list"
opts$bed = "first_panel_for_CNV.bed"
opts$outdir = "."
opts$threshold = 5

###### Check options 
######  control samples bam file 

control_bam_list <- opts$bam

control_list <- system(sprintf('cut -f7 -d"/"  %s',control_bam_list),intern = T)
#control_list <- opts$control
# reference.file <- "/data2/database/b37/human_g1k_v37.fasta"
reference.file <- "../workdir/human_g1k_v37.fasta"

cat("Step1: Reading target ranges from first.panle.bed\n")
####### read target ranges from first.panle.bed  
target.file <- opts$bed
target.df   <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end","gene"))
target      <- GRanges(seqname=target.df$seqname,IRanges(start=target.df$start+1,end=target.df$end))
mcols(target)$gene <- target.df$gene
target.sub  <- target[(target %>% start() ) %in%  (reduce(target,min.gapwidth=0) %>% start()),]

###### Filter chromosome range number less than 5   
filt_gene = which(table(target.sub$gene) >= 5)  %>% as.data.frame() %>% rownames()
target.sub  <- target.sub[target.sub$gene %in% filt_gene]
seqlevels(target.sub) <- target.sub@seqnames@values %>% as.character()

###### Filter  chromosome  with one kind of range width  
gene_list <- unique(as.character(target.sub$gene))

for (i in gene_list) {
  tmp =target.sub[target.sub$gene == i]
  if (length(tmp$gene) >= 5){
    if ( width(tmp) %>%  unique() %>% length() == 1) {
      start_tmp = tmp@ranges@start
      width_tmp = tmp@ranges@width
      width_tmp[length(width_tmp)] = width_tmp[length(width_tmp)] + 1
      tmp0 = GRanges(seqnames = tmp@seqnames,IRanges(start = start_tmp, width =  width_tmp), strand = tmp@strand)
      mcols(tmp0)$gene <- tmp$gene
      target.sub[target.sub$gene == i] = tmp0
    }else{
      target.sub[target.sub$gene == i] = tmp
    }
  }
}


cat("Step2: Read control bam files from xxx.target.realigned.rmdup.bam\n")
####### Read control bam files from xxx.target.realigned.rmdup.bam
control.list <-   control_list 
control_bam_list <- read.table(control_bam_list,header=F,stringsAsFactors = F)
control.bam <- data.frame(samples=control.list,bam=sapply(control_bam_list,function(n) sprintf("%s.target.realigned.rmdup.bam",n)),stringsAsFactors=F)
colnames(control.bam) <- c("samples","bam")


df <- control.bam
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
baseline <- counts[counts$bg > 0,]

save(baseline,file=sprintf("%s/control_%s.RData",opts$outdir,format(Sys.Date(),format="%Y%m%d")))
rm(baseline)
