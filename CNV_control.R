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
                      make_option(c("-o", "--outdir"),type="character",help="output directory")
                      )
  opts <- parse_args(OptionParser(option_list=option_list))
}


opts = list()
opts$bam = "../bam_list"
opts$bed = "first_panel_for_CNV.bed"
opts$outdir = "."
opts$threshold = 5

###### Check options 
######  control samples bam file 

control_bam_list <- opts$bam

control_list <- system(sprintf('cut -f7 -d"/"  %s',control_bam_list),intern = T)
control_list <- sapply(control_list, function(n) str_split(n[1],pattern = "\\.")[[1]][1])

#control_list <- opts$control
# reference.file <- "/home/wqh/workdir/database/human_g1k_v37.fasta"
reference.file <- "/data2/database/b37/human_g1k_v37.fasta"
cat("Step1: Reading target ranges from bed file\n")
####### read target ranges from bed file 
target.file <- opts$bed
target.df   <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end","gene"))
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

cat("Step2: Read control bam files\n")
####### Read control bam files
control.list <-   control_list 
#control.list <- read.table(control_list,header=F,stringsAsFactors = F)
control_bam_list <- read.table(control_bam_list,header=F,stringsAsFactors = F)
control.bam <- data.frame(samples=control.list,bam=control_bam_list,stringsAsFactors=F)
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
