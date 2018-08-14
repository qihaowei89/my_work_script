#! /usr/bin/R
#################################################################
### qihao
### from dbsn(hg19_avsnp147.txt) load SNPs : chr, pos , ref, aln  ,AND add column: geno
### data : 2018-5-14 
#################################################################

myf <-function(filename,output,sheet){
    library(readxl)
    library(stringr) # install.packages("stringr")
    library(tidyr) # install.packages("tidyr")
    library(parallel)
    no_cores <- detectCores()-1
    cl <- makeCluster(getOption("cl.cores", no_cores))  
    snp_list <- function(file){
      temp <- file$检测位点
      temp
    }
    file <- read_xlsx(filename,sheet = sheet)
    file <- unite(file,col="疗效或毒副作用预测(仅供参考)",c("疗效预测","毒副作用预测"),sep = ",",remove=F)
    a <- unlist(lapply(file$`疗效或毒副作用预测(仅供参考)`, function(n) str_replace(n,pattern = "NA,",replacement = ""))) 
    a <- unlist(lapply(a, function(n) str_replace(n,pattern = ",NA",replacement = ""))) 
    file$`疗效或毒副作用预测(仅供参考)` <- a 
    file <- file[c(1,2,3,4,5,8)]
    list <- snp_list(file)
    #system.time({aa = sapply(list, function(i)　system(sprintf("grep %s$ hg19_avsnp147.txt",i),intern = T))})
    #system.time({aa = parSapply(cl, list, function(i)　system(sprintf("grep %s$ hg19_avsnp147.txt",i),intern = T))})
    aa = parSapply(cl, list, function(i)　system(sprintf("grep %s$ ~/workdir/database/hg19_avsnp147.txt",i),intern = T))
    data = NULL
    for (i in names(aa)) {
      if(length(aa[[i]]) > 1){
        line1 = aa[[i]][1] %>% str_split(pattern = "\t") %>% '[['(1)
        line2 = aa[[i]][2] %>% str_split(pattern = "\t") %>% '[['(1)
        line = c(line1[6],line1[1],line1[2],line1[4],paste(line1[5],line2[5],sep = ","))
      }else{
        tmp = aa[[i]][1] %>% str_split(pattern = "\t") %>% '[['(1)
        line = c(tmp[6],tmp[1],tmp[2],tmp[4],tmp[5])
      }
      data = rbind(data,line)
    }
    fram = data.frame(file,data)
    colnames(fram) <- c("药物","基因", "检测位点","基因型","疗效或毒副作用预测(仅供参考)","等级","snp","chr","pos","ref","alt")
    dat = NULL
    for(i in 1:length(fram[,1])){
      tmp = fram[i,]
      if(is.na(tmp$alt)) {geno = NA
      }else{
        if(str_count(tmp$基因型 ,pattern =as.character(tmp$alt))==0) geno = 0
        if(str_count(tmp$基因型 ,pattern =as.character(tmp$alt))==1) geno = 0.5
        if(str_count(tmp$基因型 ,pattern =as.character(tmp$alt))==2) geno = 1
      }
      dat0 = cbind(tmp[,-7],geno)
      dat = rbind(dat,dat0)
    }
    write.table(dat,paste(output,"RS.target",sep = "."),sep = "\t",quote = F,row.names = F)
}



myf(filename = "20180509.xlsx",sheet = 1,output = "all")
myf(filename = "20180509.xlsx",sheet = 2,output = "lung")
myf(filename = "20180509.xlsx",sheet = 3,output = "jzc")
myf(filename = "1化疗药物数据库-单药简化版20180620.xlsx",sheet = 3,output = "BC") #2018-06-25 breast cancer target



