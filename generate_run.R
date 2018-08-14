
library(stringr)

args<-commandArgs(T)
# sample_file = read.table("samples_20180801.txt",stringsAsFactors = F)
sample_file = read.table(args[1],stringsAsFactors = F)
date_time = str_extract(args[1],pattern = "\\d+")

sink(sprintf("run_%s.sh",date_time))
for(i in 1:dim(sample_file)[1]){
  write.table(sprintf("snpr-gen-report -i %s" , sample_file[i,1]),quote = F,row.names = F,col.names = F)
  tmp = sprintf("snpr-gen-pdf -i  %s -d %s -n %s -C -e -u http://192.168.1.205:8083 -p CL01001  -o %s_pdf %s.gender.csv",
                sample_file[i,1],(sample_file[i,1] %>% str_split(pattern = "_",simplify = T) %>% '['(2)),sample_file[i,2],date_time,date_time)
  write.table(tmp,quote = F,row.names = F,col.names = F)
}
sink()
