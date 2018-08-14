#!/usr/bin/env Rscript
# by: qihao
# date : 2018-05-22
# for : hot spot　recovery
# 参数说明
# Options
# -f/--file       热点位点文件 　　　　
# -v/--vcf        vcf文件 　　　　　
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

library(readxl);library(tidyr);library(stringr);library(magrittr);library(optparse)

# 解析命令行
if (TRUE){
  option_list <- list(make_option(c("-f", "--file"), type="character"), 
                      make_option(c("-v", "--vcf"), type="character")
                      #make_option(c("-o","--output"), type="character")
                      )
  opts <- parse_args(OptionParser(option_list=option_list))
}

# file = read_excel("mutation_hotspot_2018_05_21.xlsx",sheet = 2)
file = read_excel(opts$file,sheet = 2)
file$`#chr` %>% str_replace(pattern="chr",replacement ="")  %>% as.numeric()-> file$`#chr`
file <- file[order(file$`#chr`,file$pos_start),]

#vcf <-read.table("Ct18050104_1011_P1080.vcf",stringsAsFactors = F)
vcf <-read.table(opts$vcf,stringsAsFactors = F)
colnames(vcf) <- c("chr","pos","a","ref","alt","b","c","d","e","f")
dat_out=NULL
for (a in 1:dim(vcf)[1]) {
  vcf_tmp <- vcf[a,][,-c(3,6,7,9,10)]
  ref = vcf_tmp$ref %>% str_split(pattern = "") %>% '[['(1)
  alt = vcf_tmp$alt %>% str_split(pattern = ",")%>% '[['(1) %>% str_split(pattern = "") 
  if (T){
    Tag =NULL
    for (i in 1:length(alt)){
      alt_tmp <- alt[[i]]
      if (!("<" %in% alt_tmp)){
        n = length(ref) - length(alt_tmp )
        if(n > 0 ){ #缺失
          INDEL_ref = ref[2:(n+1)] %>% paste0(collapse = "")
          INDEL_alt = "-"
          Tag[i] = list(c(vcf_tmp$chr,vcf_tmp$pos+1,INDEL_ref,INDEL_alt))
        }
        if(n < 0){　#插入
          INDEL_ref ="-"
          INDEL_alt = alt_tmp[2:(1-n)] %>% paste0(collapse = "")
          Tag[i] = list(c(vcf_tmp$chr,vcf_tmp$pos+1,INDEL_ref,INDEL_alt))
        }
        if(n==0){ #单碱基突变
          SNP_ref = ref
          SNP_alt = alt_tmp
          Tag[i] = list(c(vcf_tmp$chr,vcf_tmp$pos,SNP_ref,SNP_alt))
        }
      }
    }
  }
  if(!is.null(Tag)){
    for (j in 1:length(Tag)) {
      chr_t = Tag[[j]][1];pos_t = Tag[[j]][2];ref_t=Tag[[j]][3];alt_t=Tag[[j]][4]
      for(x in 1:dim(file)[1]){
        tmp = file[x,]
        if(chr_t==tmp$`#chr`& pos_t==tmp$pos_start & ref_t==tmp$ref & alt_t==tmp$alt){
          out_tmp = cbind(tmp,vcf_tmp )
          dat_out = rbind(dat_out,out_tmp[,-c(21:22)])  
        }
      }
    }
  }
}
write.table(dat_out,paste0(opts$out,".hotspot.xls"),quote = F,row.names = F)
