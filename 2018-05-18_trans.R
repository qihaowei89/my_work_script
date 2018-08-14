#!/usr/bin/env Rscript
# by: qihao
# date : 2018-05-18
# for : get intersection elements of xx.MuTect2.vcf between xx.RS.target,and reorder output file  
# 参数说明
# Options
# -v/--vcf        MuTec2.vcf文件 　　　　
# -r/--ref        参照文件 　　　　　
# -s/--sample 	  



options(warn = -1)
package_list <- c("optparse","readxl","magrittr","stringr")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#　clean enviroment object
rm(list=ls()) 

# Load essential packages
library(optparse)
library(readxl)
library(magrittr)
library(stringr)


# 解析命令行
if (TRUE){
  option_list <- list(make_option(c("-v", "--vcf"),    type="character",help = "输入xx.MuTect2.vcf"), 
                      make_option(c("-t", "--target"), type="character",help = "输入xx.RS.target文件"),
                      make_option(c("-r", "--ref"),    type="character",help = "输入排序参考文件", default = "~/workdir/handover/ctDNA_pipeline/chemical_medicine/ref_2018-05-09"),
                      make_option(c("-s", "--sample"), type="character",help = "样本结果输出路径"))
  opts <- parse_args(OptionParser(option_list=option_list,usage = "usage: %prog [options]",add_help_option = T))
}

vcf <- file(opts$vcf)
open(vcf)
sample_db =　NULL
while(T){
  line <- readLines(vcf,n = 1)
  if (length(line) == 0) break
  if(!grepl(line,pattern = "^#")){
    tmp <- str_split(line,"\t") %>% '[['(1)
    cut <- tmp[10] %>% str_split(":") %>% "[["(1) %>% "["(3) %>% as.numeric()
    if ( cut >= 0.8 ) geno = 1
    if ( cut >= 0.3 & cut < 0.8 ) geno = 0.5
    if ( cut <  0.3 ) geno = 0
    Row <- c(tmp[1],tmp[2],tmp[4],tmp[5],geno)
    sample_db <- rbind(sample_db,Row)
  }
}
close(vcf)
colnames(sample_db) <- c("chr","pos","ref","aln","geno")
sample_db  <- data.frame(sample_db,row.names = NULL,stringsAsFactors = F)
target_db  <- read.table(opts$target,header = T,stringsAsFactors = F)
CHR_sample = sample_db$chr %>% unique()
CHR_target = target_db$chr %>% unique()
db_vcf     = sapply(CHR_sample, function(i) list(sample_db[which(sample_db$chr %in% i),]))
db_target  = sapply(CHR_target, function(i) list(target_db[which(target_db$chr %in% i),]))
output = NULL
for(i in CHR_sample[CHR_sample %in% CHR_target]){
    target = db_target[[i]][which(( db_target[[i]]$pos %in% db_vcf[[i]]$pos)),]
    vcf = db_vcf[[i]][which((db_vcf[[i]]$pos) %in%  db_target[[i]]$pos ),]
    tmp = target[target$geno %in% vcf$geno,][,-(7:11)]
    if(!dim(tmp)[1]==0) output = rbind(output,tmp)
}


ref <- read.table(opts$ref,header = T, stringsAsFactors = F)[,1] %>% unique() %>% '[['(1)
out=NULL
for(i in ref){
  tmp=output[which(output$药物 %in% i),]
  if(dim(tmp)[1] > 1)
    tmp <- tmp[order(tmp$基因,tmp$检测位点),]
  out = rbind(out,tmp)
}

write.table(out,file = paste0(opts$sample,".CHEMICAL.medicine.xls"),quote = F,sep = "\t",row.names = F)




