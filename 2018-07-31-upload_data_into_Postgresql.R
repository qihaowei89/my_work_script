library(RPostgreSQL)
# install.packages("RPostgreSQL")
con = dbConnect(PostgreSQL(),host="127.0.0.1",user="postgres",password="wqh123",dbname="mydb")
dbListTables(con)
dbWriteTable(con,"mtcars",mtcars)
dbListTables(con)
dbListFields(con,"mtcars")
dbReadTable(con,"mtcars")
res = dbSendQuery(con,"select * FROM mtcars where cyl =4 order by mpg")
# You can fetch all results:
dbFetch(res)
dbClearResult(res)

# Or a chunk at a time
res <- dbSendQuery(con, "SELECT * FROM mtcars WHERE cyl = 4")
while(!dbHasCompleted(res)){
  chunk <- dbFetch(res, n = 5)
  print(nrow(chunk))
}
# Clear the result
dbClearResult(res)

# Disconnect from the database
dbDisconnect(con)


dbListTables(con)
dbListFields(con,"test0729")
dbReadTable(con,"test")
dbSendStatement(con,"drop table test0729;")
dbSendStatement(con,"create table test (id int, name varchar(20));insert into test values (1,'tom');")


conn=dbConnect(PostgreSQL(),host="192.168.1.205",port=5440,dbname="ancestry",user="postgres",password="123456")
a = dbListTables(conn)
dbDisconnect(conn)
dbReadTable(conn,a[12]) %>% tail(20)


dbSendStatement(conn,"")



library(stringr)
library(RPostgreSQL)
con = dbConnect(PostgreSQL(),host="127.0.0.1",user="postgres",password="wqh123",dbname="ancestry")
dbListTables(con)
trans_to_use = function(out_file){
  group = read.table(out_file,stringsAsFactors = F)
  group_num = str_extract(group[2:4],pattern = "\\(\\d{1}\\)")  %>% sapply(FUN = function(n) gsub(pattern = "\\((\\d)\\)", replacement = "\\1",x = n,perl = T))
  group_class = group[2:4] %>% sapply(FUN = function(n) gsub(pattern = "(.?)\\(\\d\\)", replacement = "\\1",x = n,perl = T))
  group_class[group_class == ""] = "NULL"
  group_use = data.frame(id=group[[1]],
                         final=group[[5]],
                         png=paste0(group[5],".png"),
                         group1=group_class[1],
                         num1=group_num[1],
                         group2=group_class[2],
                         num2=group_num[2],
                         group3=group_class[3],
                         num3=group_num[3],stringsAsFactors = F)
  return(group_use)
}

file = region_results_use
table = "region_results"
AddValues2Db = function(file,table){
  if(table=="region_results") {index="sample"}else{index="id"}
  if(table=="locations") break()
  # if(table=="region_results") index="sample"
  header=colnames(file)
  a = which(index %in% header)
  if(index %in% header){
    id = dbSendStatement(con,sprintf('select %s from %s;',index,table)) %>%  dbFetch(n = -1)  %>% unique() %>% '[['(1)
    if(any(id %in% unique(file[,a]))){
      if(any(grepl("_",id))){
        id_tmp=unique(file[,a]) %>% str_split(pattern = "_",simplify = T) %>% '['(2)
        dbSendStatement(con,sprintf('delete from %s where %s like %s%s%s ',table,index,"'%",id_tmp,"'"))
      }else{
        dbSendStatement(con,sprintf('delete from %s where %s = %s ',table,index,unique(file[,a])))
      }
    }
  }
  dbWriteTable(con,table,file,append=T,row.names=F)
}
if(F){
    dbSendStatement(con,
         "CREATE TABLE REGION_RESULTS( ID SERIAL PRIMARY KEY NOT NULL,
                    SAMPLE  TEXT NOT NULL,
                    REGION  TEXT NOT NULL,
                    PERCENT INTEGER);
          CREATE TABLE LOCATIONS( ETHNIC  TEXT  PRIMARY KEY   NOT NULL,
                    LOCAL   TEXT NOT NULL,
                    OCEANIA TEXT NOT NULL);
          CREATE TABLE NAMES ( LOCATION TEXT PRIMARY KEY   NOT NULL,
                    CHINESE TEXT NOT NULL);
          CREATE TABLE CONTENTS( ID INTEGER PRIMARY KEY   NOT NULL,
                    COVER_PAGE TEXT NOT NULL,
                    REGION_PAGE TEXT NOT NULL,
                    OUTLINE TEXT NOT NULL);
          CREATE TABLE Y_DESCRIPTION( GROUPS TEXT PRIMARY KEY NOT NULL,
                    DESCRIPTION TEXT NOT NULL);
          CREATE TABLE SAMPLE_GENDER( SAMPLE TEXT PRIMARY KEY NOT NULL,
                    GENDER TEXT NOT NULL);
          CREATE TABLE Y_GROUP (ID TEXT PRIMARY KEY     NOT NULL,
                    FINAL          TEXT     NOT NULL,
                    PNG            TEXT     NOT NULL,
                    group1         TEXT    ,
                    num1            INT     NOT NULL,
                    group2         TEXT    ,
                    num2            INT     NOT NULL,
                    group3         TEXT    ,
                    num3            INT     NOT NULL
                    );
          CREATE TABLE MT_GROUP(ID TEXT PRIMARY KEY     NOT NULL,
                    FINAL          TEXT     NOT NULL,
                    PNG            TEXT     NOT NULL,
                    group1         TEXT    ,
                    num1            INT     NOT NULL,
                    group2         TEXT    ,
                    num2            INT     NOT NULL,
                    group3         TEXT    ,
                    num3            INT     NOT NULL);")
}

dbListTables(con)
region_results = read.table("/home/wqh/B10_10003858.region.out",sep = "\t",header = T,stringsAsFactors = F)
region_num = round(as.numeric(region_results[,2][-58])*100)
region_results_use = data.frame(sample="B10_10003858",region=as.character(region_results[-58,1]),percent=region_num,stringsAsFactors = F)

content = read.table("/home/wqh/content_2.csv",sep = ",",header = T,stringsAsFactors = F)


names_db = read.table("/home/wqh/ethnic_name.csv",sep=",",header=F,stringsAsFactors = F)
head(names_db)
colnames(names_db) = c("ethnic","ethnic_ch","local","local_ch","oceania","oceania_ch")
names_db$ethnic %<>% sapply(FUN=function(n) str_replace_all(string = n,pattern = " ",replacement = "_"),simplify = T)
locations = names_db[,c(1,3,5)]
names = names_db[,c(1,2)] 
colnames(names) = c("location","chinese")

y_group  = trans_to_use("/home/wqh/B10_10003858.Y.out")
MT_group = trans_to_use("/home/wqh/B10_10003858.MT.out")




dbReadTable(con,"region_results")
dbWriteTable(con,"region_results",region_results_use,append=T,row.names=F)
AddValues2Db(file = region_results_use,table = "region_results")
dbReadTable(con,"region_results")

# dbSendStatement(con,"drop table y_group")
dbWriteTable(con,"content",content ,append=T,row.names=F)
dbReadTable(con,"content")

AddValues2Db(file =content, table="content")

dbWriteTable(con,"locations",locations,append=T,row.names=F)
AddValues2Db(file =locations, table="locations")
dbReadTable(con,"locations")

dbWriteTable(con,"names",names,append=T,row.names=F)
dbReadTable(con,"names")

AddValues2Db(file =y_group, table="y_group")
dbWriteTable(con,"y_group",y_group,append=T,row.names=F)
dbReadTable(con,"y_group")


dbWriteTable(con,"mt_group",MT_group,append=T,row.names=F)
dbReadTable(con,"mt_group")



