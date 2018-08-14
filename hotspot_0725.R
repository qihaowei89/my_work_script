#!/usr/bin/env Rscript
# by: qihao
# date : 2018-05-22,2018-07-09_1,2018-07-09_2,2018-07-25
# for : hot spot　recovery
# 参数说明
# Options
# -d/--dir        热点位点文件所在文件夹 　　　　
# -v/--vcf        vcf文件 　　　　　
# -o/--outdir     输出文件
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


if(TRUE){
  option_list <- list(make_option(c("-d", "--dir"), type="character",help = "mutation_hotspot Excel file dir"), 
                      make_option(c("-v", "--vcf"), type="character"), 
                      make_option(c("-o","--outdir"), type="character"))
  opts <- parse_args(OptionParser(option_list=option_list))
}

file = read_excel(sprintf("%s/mutation_hotspot.xlsx",opts$dir),sheet = 2)


file$`#chr` %>% str_replace(pattern="chr",replacement ="")  %>% as.numeric() -> file$`#chr`
file <- file[order(file$`#chr`,file$pos_start),]

vcf <-read.table(opts$vcf,stringsAsFactors = F)
colnames(vcf) <- c("chr","pos","a","ref","alt","b","c","d","e","f")

dat_out=NULL
for (a in 1:dim(vcf)[1]) {
  vcf_tmp <- vcf[a,][,-c(3,6,7,9,10)]
  ref = vcf_tmp$ref 
  alt = vcf_tmp$alt %>% str_split(pattern = ",")%>% '[['(1) 
  alt = alt[!(alt %in% c("<*>"))] 
  Tag =NULL
  if (!length(alt)==0) {
    for (i in 1:length(alt)){
      alt_tmp <- alt[i]
      if (!("<" %in% alt_tmp)){
        n = nchar(ref) - nchar(alt_tmp)
        if(n >  0 ){ #缺失
          ref_tmp = ref %>% str_split(pattern = "") %>% '[['(1)
          INDEL_ref = ref_tmp[(nchar(alt_tmp)+1):nchar(ref)] %>% paste0(collapse = "")
          INDEL_alt = "-"
          Tag[i] = list(c(vcf_tmp$chr,vcf_tmp$pos+1,INDEL_ref,INDEL_alt))
        }
        if(n <  0){ #插入
          alt_tmp_tmp=alt_tmp %>% str_split(pattern = "") %>% unlist()
          INDEL_ref ="-"
          INDEL_alt = alt_tmp_tmp[(nchar(ref)+1):nchar(alt_tmp)] %>% paste0(collapse = "")
          Tag[i] = list(c(vcf_tmp$chr,vcf_tmp$pos+1,INDEL_ref,INDEL_alt))
        }
        if(n == 0){ #单碱基突变
          SNP_ref = ref
          SNP_alt = alt_tmp
          Tag[i] = list(c(vcf_tmp$chr,vcf_tmp$pos,SNP_ref,SNP_alt))
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
}

split = opts$vcf %>% str_split(pattern = "/",simplify = T)
filename = split[length(split)] %>% str_split(pattern = "\\.",simplify = T) %>%'['(1)
write.table(dat_out,paste0(opts$outdir,"/",filename,".hotspot.xls"),sep = "\t",quote = F,row.names = F)
