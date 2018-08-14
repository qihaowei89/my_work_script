#!/usr/bin/env Rscript
# by: qihao
# date : 2018-05-18
# for : reorder output file 
# 参数说明
# Options
# -f/--file       待排序文件 　　　　
# -r/--ref        参照文件 　　　　　

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
  option_list <- list(make_option(c("-f", "--file"), type="character"), make_option(c("-v", "--vcf"),  type="character"))
  opts <- parse_args(OptionParser(option_list=option_list))
}
#file = read_excel("mutation hotspot.xlsx",sheet = 2)
file = read_excel(opts$file,sheet = 2)
file$`#chr` %>% str_replace(pattern="chr",replacement ="")  %>% as.numeric()-> file$`#chr`
file <- file[order(file$`#chr`,file$pos_start),]
# pos = file[,c(1,2,3)]
# colnames(pos) = c("chr","start","end") 
# write.table(pos,file = "ref_sorted",quote = F,row.names = F,col.names = F)
vcf = read.table(opts$vcf,stringsAsFactors = F)[,1:8][,c(1,2,4,5,8)]
#vcf = read.table("Ct18050110_1011_T1084.vcf",stringsAsFactors = F)[,1:8][,c(1,2,4,5,8)]
colnames(vcf) = c("chr","pos","ref","aln","meta")
chr = intersect(unique(vcf$chr),unique(file$`#chr`))
vcf_db = vcf[vcf$chr %in% chr,]
file_db = file[file$`#chr` %in% chr,]
dat_db=NULL
###
for(i in seq(dim(file_db)[1])){
  file_tmp = file_db[i,]
  ###是否为插入或者缺失
  if(file_tmp$ref == "-"|file_tmp$alt == "-"){
    if(file_tmp$alt == "-"){#缺失
      for (j in seq(dim(vcf_db)[1])){
        vcf_tmp = vcf_db[j,]
        REF = vcf_tmp$aln %>% str_split(pattern = ",") %>% '[['(1)
        vcf_pos = vcf_tmp$pos+1
        if(vcf_tmp$ref %>% str_split("") %>% '[['(1) %>% '['(1) == REF) REF = "-"
        vcf_ref = vcf_tmp$ref %>% str_split("") %>% '[['(1) %>% '['(-1) %>% paste0(collapse = "")
        if(vcf_tmp$chr==file_tmp$`#chr` & vcf_pos==file_tmp$pos_start & vcf_ref==file_tmp$ref & (file_tmp$alt %in% REF)){
          tmp = cbind(file_tmp,(vcf_tmp %>% separate(col = "meta",into = letters[1:13],sep = ";") %>% '['(1:8)))
          dat_db = rbind(dat_db,tmp)
        }

      }
    }
    if(file_tmp$ref == "-"){#插入
      # for (j in seq(dim(vcf_db)[1])){
      #   vcf_tmp = vcf_db[j,]
      #   REF = vcf_tmp$aln %>% str_split(pattern = ",") %>% '[['(1)
      #   vcf_pos = vcf_tmp$pos+1
      #   if(length(vcf_tmp$ref %>% str_split("") %>% '[['(1) ) > 1 ){
      #     if(vcf_tmp$ref %>% str_split("") %>% '[['(1) %>% '['(1) == REF) REF = "-"
      #     vcf_ref = vcf_tmp$ref %>% str_split("") %>% '[['(1) %>% '['(-1) %>% paste0(collapse = "")
      #   }
      #   if(vcf_tmp$chr==file_tmp$`#chr` & (vcf_tmp$pos+1)==file_tmp$pos_start & vcf_ref==file_tmp$ref & (file_tmp$alt %in% REF)){
      #     tmp = cbind(file_tmp,(vcf_tmp %>% separate(col = "meta",into = letters[1:13],sep = ";") %>% '['(1:8)))
      #     dat_db = rbind(dat_db,tmp)
      #   }
      # }
      
    }
  }else{##为普通单碱基突变
    for (j in seq(dim(vcf_db)[1])){
      vcf_tmp = vcf_db[j,]
      REF = vcf_tmp$aln %>% str_split(pattern = ",") %>% '[['(1)
      if(vcf_tmp$chr==file_tmp$`#chr`& vcf_tmp$pos==file_tmp$pos_start&vcf_tmp$ref==file_tmp$ref&(file_tmp$alt %in% REF)){
        tmp = cbind(file_tmp,(vcf_tmp %>% separate(col = "meta",into = letters[1:13],sep = ";") %>% '['(1:8)))
        dat_db = rbind(dat_db,tmp)
      }
    }
  }
}

write.table(dat_db,quote = F,col.names = T,row.names = F)

  


