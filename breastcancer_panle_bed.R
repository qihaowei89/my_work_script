library(stringr)
library(tidyr)

file0 = read.table("../wqh/下载/designed-probe-coords.bed",stringsAsFactors = F)
file1 = read.table("../wqh/下载/designed-probe-coords[1].bed",stringsAsFactors = F)

file = rbind(file0,file1)


head(file)
file = separate(file,col="V4",into =c("a","b","c","d") ,sep = "_")[,c(-4,-5,-7,-8,-9)]

file$gene = str_extract(file$c,pattern = "[A-Z0-9]+") 



file[439:440,]$gene <- "TERT"
file[581,]$gene <- "CDKN2A"
file[578:583,]

head(file[,-4])

a= file

file = file[,-4]
file$V1 = str_extract(file$V1,pattern = "[0-9A-Z]+") %>% as.numeric() 

file$V1[which(is.na(file$V1))] <- 23

file = file[order(file$V1),]
file
file$V1[which(file$V1 == 23)] = "X"
tail(file)

file = read.table("first_panel_for_CNV.bed",stringsAsFactors = F)
write.table(file,file = "first_panel_for_CNV.bed",sep = "\t", quote = F,col.names = F,row.names = F)
getwd()



?GRanges

gr0 <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               IRanges(1:10, width=10:1))
gr0

names(gr0) <- head(letters, 10)

strand(gr0) <- Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2))
mcols(gr0)$score <- 1:10
mcols(gr0)$GC <- seq(1, 0, length=10)
gr0
