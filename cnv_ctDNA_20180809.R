#!/usr/bin/env Rscript
# by: qihao
# date : 2018-06-06, 2018-06-27,2018-07-16

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


if (TRUE){
  option_list <- list(make_option(c("-B", "--bam"),type="character",help="input bam file"), 
                      make_option(c("-b", "--bed"),type="character",help="target region bed file"),
                      make_option(c("-s", "--sample"),type="character",help="sample name"),
                      make_option(c("-c", "--controlset"),type="character",help="control set"),
                      make_option(c("-o", "--outdir"),type="character",help="output directory"),
                      make_option(c("-t", "--threshold"),type="numeric",default=3,help="copy number threshold to output [default %default]"))
  opts <- parse_args(OptionParser(option_list=option_list))
}


opts= list()
# opts$bam = "/home/wqh/workdir/bam/Ct1807_1113_20020702.target.bam"
# opts$bam = "/home/wqh/Ct1807_1113_10020395.target.rmdup.realigned.bam"
# opts$bam = "/home/wqh/Ct1807_1113_20020858.target.rmdup.realigned.bam"
opts$bam = "/home/wqh/workdir/bam/Ct1807_1113_CNV_1.target.bam"

opts$bed = "/home/wqh/ctDNA_103genes_20180702.bed"
opts$sample = "Ct1807_1113_CNV_1"
# opts$sample = "Ct1807_1113_10020395"
# opts$sample = "Ct1807_1113_CNV_1"
# opts$sample = "Ct1807_1113_20020858"

opts$controlset = "/home/wqh/workdir/bam/cnv_controlset_ctDNA.RData"
opts$outdir = "."
opts$threshold = 3



###### Check options 
if(length(opts) < 7){
  cat ("Usage: Rscript cnv_ctDNA.R -B <Input bamFile> -b <Input bedFile> -s <Input samples list> -c <Control set> -o <Output directory> -t <Threshold>\n")
  cat ("Example: Rscript cnv_ctDNA.R -B ct_1013_HD786_1.target.rmdup.realigned.bam -b first.panel.bed -s ct_1013_HD786_1 -c cnv_controlset_ctDNA.RData -o ./ -t 3\n")
  cat ("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat ("Version: v1.0\n")
  cat ("Date: 2018-06-28\n")
}

cat("Step1: Reading target ranges from bed file\n")
####### read target ranges from bed file  
target.file <- opts$bed
target.df   <- read.delim(target.file,header=FALSE,col.names=c("seqname","start","end","gene"),stringsAsFactors = F)
target      <- GRanges(seqname=target.df$seqname,IRanges(start=target.df$start+1,end=target.df$end))
target.sub  <-target
# target.sub  <- reduce(target,min.gapwidth=0)
target.gene <- target.sub
mcols(target)$gene <- target.df$gene 

target@ranges %>% width() %>% table

# subdivideGRanges()

#gene list 
target.gene <- target[(target %>% start()) %in% (reduce(target,min.gapwidth=0) %>% start()),]
target.gene <- target
target.gene@ranges %>% width()

###### Filter chromosome range number less than 5   
# target.sub  <- target.sub[target.sub@seqnames %in% target.sub@seqnames@values[target.sub@seqnames@lengths >= 5]]
# target.gene <- target.gene[target.sub@seqnames %in% target.sub@seqnames@values[target.sub@seqnames@lengths >= 5]]
# # seqlevels(target.sub) <- target.sub@seqnames@values %>% as.character()

###### Filter  chromosome  with one kind of range width  
# seqs_list <- target.sub@seqnames@values
# for (i in seqs_list) {
#   tmp =target.sub[seqnames(target.sub) == i]
#   if (tmp@seqnames@lengths >= 5){
#     if ( width(tmp) %>%  unique() %>% length() == 1) {
#       start_tmp = tmp@ranges@start
#       width_tmp = tmp@ranges@width
#       width_tmp[length(width_tmp)] = width_tmp[length(width_tmp)] + 1
#       tmp = GRanges(seqnames = i,IRanges(start = start_tmp, width =  width_tmp), strand = tmp@strand)
#       target.sub[seqnames(target.sub) == i] = tmp
#     }else{
#       target.sub[seqnames(target.sub) == i] = tmp
#     }
#   }
# }

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
counts$gene <- target.gene$gene
counts$log.bg <- baseline$log.bg
counts$GC <- baseline$GC 
counts$GC.sq <- baseline$GC.sq
counts$width <- baseline$width

counts[which(counts@seqnames==17 & counts$gene =="ERBB2")] 

counts$bg <- baseline$bg
base_data = baseline@elementMetadata[,c(1,2,3,17)] 

cor(base_data$ct_1013_bai_1,base_data$bg)
cor(base_data$ct_1013_bai_2,base_data$bg)

rm(baseline)
fit = ExomeCopy(counts,sample.name = "Ct1807_1113_20020702")
counts$copy_num = round(((fit@O.norm / width(counts@ranges)) * 120) * 2 ,digits = 1)

# fit@O.norm[which(counts@seqnames==17 & counts$gene =="ERBB2")]*2  %>% mean()
# counts$copy_num[which(counts@seqnames==17)][12:28] %>% mean()

# counts$raw_reads_bg <-round( (counts@elementMetadata[[1]] / counts@elementMetadata[[6]]) / median(counts@elementMetadata[[1]] / counts@elementMetadata[[6]]) * 2, digits = 1)
M = (counts@elementMetadata[[1]] / counts$bg) 


copy_num_raw_new = round((M / median(M) )[which(counts@seqnames==17)][22:39] * 2 ,digits = 1) %>% mean()
  
counts$raw_reads_bg <-round( ( counts@elementMetadata[[1]]/ (summary(counts@elementMetadata[[1]])[3:4] %>% mean())) * 2, digits = 1)

counts$raw_reads_bg <- round((counts$raw_reads_bg / width(counts@ranges)) *120,digits = 1)

width(counts@ranges)

#save(counts,file="sample.RData")

###### Plotting CNVs into pdf, Compiling and filtering results, Writeing CNVs results into files
system(sprintf("mkdir -p %s/%s.CNVs.plots",opts$outdir,opts$sample))
cat("Step4: Plotting CNVs into pdf\n       Writing CNVs results into files\n")

genes <- counts$gene %>% unique()

out = NULL
out_raw = NULL
for (gene in genes) {
  tmp = counts[counts$gene == gene] 
  tmp$width = tmp@ranges  %>% width()
  
  x = apply(ranges(tmp) %>% as.matrix(), 1, function(n) {n[1]+n[2]/2})
  y = tmp$copy_num
  MAX = max(y)
  copy_num = tmp$copy_num %>% mean() %>% round(digits = 1)
  copy_num = sum(tmp$copy_num * tmp$width) /sum(tmp$width)
  ((tmp$copy_num / tmp$width) * 120 ) %>% mean() %>% round(digits = 1)
  
  copy_num_raw = ((tmp$raw_reads_bg / tmp$width)*120)  %>% mean() %>% round(digits = 1)
  
  copy_num_raw_new = 
  (sum(tmp@elementMetadata$raw_reads_bg[[1]] * tmp@elementMetadata$width[[1]]) / sum(tmp@elementMetadata$width[[1]])) * 2
  
  if(copy_num >= opts$threshold){
    out_tmp = DataFrame(chromosome = tmp@seqnames@values,gene = gene,ranges = IRanges(start = start(tmp)[1],end=end(tmp)[length(tmp)]), copy.num = copy_num,nranges = length(tmp))
    out = rbind(out,out_tmp)
    pdf(sprintf("%s/%s.CNVs.plots/%s_%s.pdf",opts$outdir,opts$sample,opts$sample,gene),width = 10,height = 6)
    plot(x,y ,xlab = sprintf("Position of chromosome %s",tmp@seqnames@values),ylab="Copy number", main = sprintf("Copy number of gene  %s",gene),ylim=c(0,MAX+2))
    abline(h = 2, col ="red",lty = 3,type = "a")
    abline(h = copy_num, col ="blue",lty = 4,type = "b")
    dev.off()
  }
  
  if(copy_num_raw >= opts$threshold){
    out_raw_tmp = DataFrame(chromosome = tmp@seqnames@values,gene = gene,ranges = IRanges(start = start(tmp)[1],end=end(tmp)[length(tmp)]), copy.num = copy_num_raw,nranges = length(tmp))
    out_raw = rbind(out_raw,out_raw_tmp)
  }
}


if(length(out) != 0){
  #reorder output file 
  Chromosome=as.character(out$chromosome) %>% as.numeric() 
  Chromosome[is.na(Chromosome)] = 23
  Out = DataFrame(Sample = opts$sample,Chromosome=Chromosome,Start=start(out$ranges),End=end(out$ranges),Gene=out$gene,Segment_length=width(out$ranges),Copy_number=out$copy.num,Bin_count=out$nranges)
  Out = Out[order(Out$Chromosome,Out$Start),]
  Out$Chromosome[Out$Chromosome == 23] ="X"
  write.table(Out,sprintf("%s/%s.CNV.xls",opts$outdir,opts$sample),quote = F,sep = "\t",row.names = F)
}else{
  system(sprintf("touch %s/%s.CNV.xls",opts$outdir,opts$sample))
}


if(length(out_raw) != 0){
  #reorder output file 
  Chromosome=as.character(out_raw$chromosome) %>% as.numeric() 
  Chromosome[is.na(Chromosome)] = 23
  Out_Raw = DataFrame(Sample = opts$sample,Chromosome=Chromosome,Start=start(out_raw$ranges),End=end(out_raw$ranges),Gene=out_raw$gene,Segment_length=width(out_raw$ranges),Copy_number=out_raw$copy.num,Bin_count=out_raw$nranges)
  Out_Raw = Out_Raw[order(Out_Raw$Chromosome,Out_Raw$Start),]
  Out_Raw$Chromosome[Out_Raw$Chromosome == 23] ="X"
  write.table(Out_Raw,sprintf("%s/%s.CNV_raw.xls",opts$outdir,opts$sample),quote = F,sep = "\t",row.names = F)
}else{
  system(sprintf("touch %s/%s.CNV_raw.xls",opts$outdir,opts$sample))
}


### Exomecopy  calculation cnv -- all results 
all_out = DataFrame(sample = opts$sample,chromosome =counts@seqnames,ranges = counts@ranges,gene = counts$gene, copy.num = counts$copy_num,copy.num_m = counts$raw_reads_bg)
write.table(all_out,sprintf("%s/%s.copy_num.txt",opts$outdir,opts$sample),quote = F,sep = "\t",row.names = F)
# write.table(all_out_raw,sprintf("%s/%s.copy_num_raw.txt",opts$outdir,opts$sample),quote = F,sep = "\t",row.names = F)

rm(list = ls())

cat("Complete!")
