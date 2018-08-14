#!/usr/bin/env Rscript
# by: qihao
# date : 2018-05-28
# for :  filte somtic.xls file
# 参数说明
# Options
# -f/--file       热点位点文件 　　　　
# -o/--output     输出文件
options(warn = -1)
package_list <- c("optparse","readxl","magrittr","stringr","tidyr")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# clean enviroment object
rm(list=ls()) 
library(readxl);library(magrittr);library(optparse)

# 解析命令行
if (TRUE){
  option_list <- list(make_option(c("-f", "--file"), type="character"),
                      make_option(c("-s", "--sites"), type="character"))
  opts <- parse_args(OptionParser(option_list=option_list))
}

file = read.table(opts$file,header = F,sep = "\t",stringsAsFactors = F)
filters = c("intronic","UTR5","UTR3","ncRNA_intronic","ncRNA_exonic","intergenic")
str = "ABL1 /AKT1 /ALK /APC /AR /ARAF /ATM /BRAF /CCND1 /CDH1 /CDK4 /CDK6 /CDKN1A /CDKN2A /CTNNB1 /DDR2 /EGFR /ERBB2 /ERBB3 /ERBB4 /ESR1 /FBXW7 /FGFR1 /FGFR2 /FGFR3 /FGFR4 /FLT3 /GNA11 /GNAQ /GNAS /HRAS /IDH1 /IDH2 /JAK1 /JAK2 /JAK3 /KDR /KIT /KRAS /MAP2K1 /MAP2K2 /MAPK1 /MTOR /NF1 /NF2 /NRAS /NTRK1 /NTRK2 /NTRK3 /PDGFRA /PDGFRB /PIK3CA /POLE /PTCH1 /PTEN /RB1 /RET /ROS1 /SMAD4 /SMARCA4 /SMO /STAT3 /STK11  /TP53 /TSC1 /TSC2 /VHL"
gene_set =  sapply(str_split(str,pattern = "/"), function(n) trimws(n),simplify = F) %>% '[['(1)
file_fiter = file[(file$V11 %in% gene_set)& !(file$V12 %in% filters) |( file$V11 == "MET" & !(file$V12 %in% filters[-1]))|( file$V11 == "TERT" & !(file$V12 %in% filters[c(-2,-3)])), ]
dele_site = read_xlsx(opts$sites,sheet = 1,col_names = F)                  
out = file_fiter[!(file_fiter$V1 %in% dele_site$X__1 & file_fiter$V2 %in% dele_site$X__2 & file_fiter$V3 %in% dele_site$X__3& file_fiter$V4 %in% dele_site$X__4 &file_fiter$V5 %in% dele_site$X__5),]
write.table(out,opts$file,quote = F,sep = "\t",row.names = F,col.names = F)
