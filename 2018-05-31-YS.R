

library(magrittr)
read.csv("../YS_20180531.csv")[,1] %>% unique() %>%  write.table("../YS_20180518.list",quote = F,row.names = F,col.names = F,sep = "\t")
install.packages("readxl")

library(readxl)
library(stringr) #install.packages("stringr")
setwd("Desktop/pdf/")
a = read_xlsx("~/Downloads/叶酸结果.xlsx",sheet = 1)
name = a[,1] %>% na.omit() %>% unique()
ID = a[,2] %>% unique()
c = cbind(name,ID)
filename_1 = list.files()
filename = filename_1
for(a in 1:length(c[,2])){
  i = c[a,2]
  string = filename[grep(pattern = i,filename)] %>% str_split(pattern = "\\.") %>% '[['(1)
  filename[grep(pattern = i,filename)] = sprintf("%s-%s.pdf",string[1],c[a,1])
}
file.rename(from = filename_1, to = filename)
getwd()






