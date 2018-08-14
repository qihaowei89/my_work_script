# library(RCurl)
library(stringr)
library(xml2)
library(rvest) 
library(magrittr)
library(parallel)

source =read_html("/home/wqh/桌面/Clinical Trials Search - My Cancer Genome.html")
a = html_nodes(source,"body div.main-content #clinical-trials-index p")  %>% html_text() %>% grep(pattern = "NCT ID:") 
b = html_nodes(source,"body div.main-content #clinical-trials-index p")  %>% html_text() %>% grep(pattern = "biomarkers:") 
NCT_list = html_nodes(source,"body div.main-content #clinical-trials-index p")  %>% html_text() %>% '['(a) %>% gsub(pattern = "NCT ID: ",replacement = "")
Biomakers = html_nodes(source,"body div.main-content #clinical-trials-index p")  %>% html_text(trim = T) %>% '['(b) %>% gsub(pattern = "biomarkers:",replacement = "") 
pp = data.frame(NCT_ID=NCT_list,biomakers=Biomakers,stringsAsFactors = F)



url = "https://beta.mycancergenome.org/content/clinical_trials"
site = paste(url,NCT_list ,sep = "/")

get_items = function(site){
  library(RCurl)
  library(stringr)
  library(REIDS)
  library(xml2)
  library(rvest) 
  library(magrittr)
  library(parallel)
  webpage <- read_html(site[i])
  cell = html_nodes(webpage,css = ".clinical-trial-detail .row p" ) 
  Drug = html_nodes(webpage,css = "table tr td") %>% '['(1) %>%html_text(trim = T)
  Phase = cell[3]  %>% html_text() %>% trim
  NCT_ID = html_nodes(webpage,css = ".clinical-trial-detail .body h4") %>% '['(1) %>% html_text(trim = T)
  Recruiting_Status =cell[2]  %>% html_text(trim = T) 
  Related_Conditions=html_nodes(webpage,css = ".clinical-trial-detail .row ul .condition") %>% html_text(trim = T)  %>% paste(collapse = ";")
  # Biomarkers
  Description = cell[1] %>% html_text()%>% str_split(pattern = "\n",simplify = T) %>%trim  %>% paste(collapse = "")
  URL=cell[4] %>% html_text(trim = T)
  if(length(Drug) == 0) Drug = "N/A"
  data_out = data.frame(Drug=Drug,Phase =Phase,NCT_ID =NCT_ID,Recruiting_Status=Recruiting_Status,Related_Conditions_or_Diseases=Related_Conditions,Biomakers=1,Description=Description,URL=URL,stringsAsFactors = F)
  return(data_out)
} 


# no_cores <- detectCores()-1
# cl <- makeCluster(getOption("cl.cores", no_cores)) 
# drug_info = parSapply(cl, site[1:1000], get_items,simplify = T)
# stopCluster(cl)
# out = as.data.frame(t(drug_info),stringsAsFactors = F)

write.table(cbind("Drug","Phase","NCT_ID","Recruiting_Status","Related_Conditions_or_Diseases","Biomakers","Description","URL"),"My_cancer_genome_db.xls",quote = F,append = T,col.names = F,sep = '\t',row.names = F)
for (i in 1:846) {
  print(i)
  Biomakers_tmp = pp[i,]
  items_tmp = get_items(site[i])
  if (Biomakers_tmp$NCT_ID==items_tmp$NCT_ID) {
    items_tmp$Biomakers = Biomakers_tmp$biomakers
    write.table(items_tmp,"My_cancer_genome_db.xls",quote = F,append = T,col.names = F,sep = '\t',row.names = F)
  }
}




