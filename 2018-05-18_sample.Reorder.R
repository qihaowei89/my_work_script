#!/usr/bin/env Rscript
# by: qihao
# date : 2018-05-18
# for : reorder output file 
# 参数说明
# Options
# -f/--file       待排序文件 　　　　
# -r/--ref        参照文件 　　　　　

options(warn = -1)
package_list <- c("optparse","readxl","magrittr")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# clean enviroment object
rm(list=ls()) 

#  Load essential packages
library(optparse)
library(readxl)
library(magrittr)

# 解析命令行
if (TRUE){
  option_list <- list(make_option(c("-f", "--file"), type="character"), make_option(c("-r", "--ref"),  type="character"))
  opts <- parse_args(OptionParser(option_list=option_list))
}
ref <- read.table(file =opts$ref,header = T,stringsAsFactors = F)[,1] %>% unique() 
# ref <- read_xlsx(opts$ref,sheet = 1)[,1] %>% unique() %>% '[['(1)
#file <- read_xlsx(opts$file,sheet = 1)
file <- read.table(opts$file,header = T,stringsAsFactors = F)
out=NULL
for(i in ref){
  tmp=file[which(file$药物 %in% i),]
  if(dim(tmp)[1] > 1)
    tmp <- tmp[order(tmp$基因,tmp$检测位点),]
  out = rbind(out,tmp)
}
write.table(out,file = opts$file,quote = F,sep = "\t",row.names = F)


