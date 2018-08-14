library(RMySQL)
library(magrittr)
##--------------Step 1 Connecte to  MySQL Database--------------
# con=dbConnect(dbDriver("MySQL"), dbname = "rmysql", user="root",password="123", host="127.0.0.1", port=3306,client.flag= CLIENT_MULTI_STATEMENTS)
con=dbConnect(RMySQL::MySQL(),dbname="rmysql",user="root",password="123",host="127.0.0.1",port=3306,client.flag=CLIENT_MULTI_STATEMENTS)
dbGetInfo(con) # show database infomations


dbReadTable(con,"t_blog2")
dbRemoveTable(con,"t_blog2")

##--------------Step 2 Creat a new empty table e.g:t_blog--------------
dbSendStatement(con,"CREATE TABLE t_blog2( 
                id INT PRIMARY KEY AUTO_INCREMENT, 
                title varchar(12) NOT NULL, 
                author varchar(12) NOT NULL,  
                length int NOT NULL, 
                create_date timestamp NOT NULL DEFAULT now()
                )ENGINE=INNODB DEFAULT CHARSET=utf8;")

dbListTables(con)  # show all tables  in Databases 
dbListFields(con,name ="t_blog2") # show one table's colnames

##--------------Add content into table t_blog2--------------
dbSendQuery(con,"SET NAMES UTF8")
# first approach use "dbSendStatement"
dbSendStatement(con,"INSERT INTO t_blog3(title,author,length) values('R插入的新文章','Conan',50)")
# second approach use "dbWirteTable" append=T
b = data.frame(title="R插入的新文章",author="Conan",length=50,stringsAsFactors = F)

dbWriteTable(con,"t_blog2",b,append=T,row.names=F)

dbReadTable(con,"t_blog2")
Encoding(b$title)



dbSendQuery(con,'SET NAMES UTF8')
dbReadTable(con,"t_blog")
dbSendQuery(con,"INSERT INTO t_blog(title,author,length) values('R插入的新文章23','Conan',50)")
# dbSendStatement(con,"INSERT INTO t_blog(title,author,length) values('R插入的新文章2','Conan',50)")


# dbSendStatement(con, "SELECT * FROM t_blog") %>% dbFetch(n=-1)
dbReadTable(con,"t_blog2")
b = a[1,c(2,3,4)]
b$title = "hello world~"
b
dbWriteTable(con,"t_blog2",b,append=T,row.names=F,encode="UTF8")
dbSendQuery(con,"SET NAMES UTF8")
dbReadTable(con,"t_blog2")
dbRemoveTable(con,"t_blog2")
dbDisconnect(con)

?dbWriteTable





library(RMySQL)
con <- dbConnect(dbDriver("MySQL"), dbname = "rmysql", user="root",password="123", host="127.0.0.1", port=3306)
dbListTables(con)
#写数据库表
fruits <-data.frame(id=1:5,name=c("苹果","香蕉","梨子","玉米","西瓜"),price=c(8.8,4.98,7.8,6,2.1),status=c("无","打折","无","售罄","批发"),s=date())

dbWriteTable(con,"fruits",fruits)
dbReadTable(con,"fruits")
dbRemoveTable(con,"fruits")


#写数据表，覆盖追加
testA <-data.frame(id=1:6,e=c("a","b","c","d","e","f"),c=c("我","的","世","界","变","得"))
testB <-data.frame(id=7:13,e=c("g","h","i","j","k","l","m"),c=c("奇","妙","跟","难","以","言","喻"))
#直接写testA写入test表中
dbWriteTable(conn,"test",testA,row.names=F)
dbReadTable(conn,"test")
#追加写testB追加在test表后
dbWriteTable(conn,"test",testB,append=T,row.names=F)
dbReadTable(conn,"test")
#覆盖写testB覆盖test表
dbWriteTable(conn,"test",testB,overwrite=T,row.names=F)
dbReadTable(conn,"test")



#用SQL语句查询dbGetQuery()和dbSendQuery()两种方法
dbGetQuery(conn, "SELECT * FROM fruits limit 2")

res <- dbSendQuery(conn, "SELECT *FROM fruits")
data <- dbFetch(res, n=2) #取前2条数据，n=-1时是获取所有数据
data
data <- dbFetch(res, n=-1) #取余下所有数据

dbClearResult(res)
dbDisconnect(conn) #断开连接


#用SQL语句批量查询
con1 <- dbConnect(MySQL(),host="host",dbname="mysql",user="root",password="123",port = 3306,client.flag= CLIENT_MULTI_STATEMENTS) #client.flag设置这样支持批量查
con <-dbConnect(dbDriver("MySQL"), dbname = "mysql", user="root",password="123", host="127.0.0.1", port=3306,client.flag= CLIENT_MULTI_STATEMENTS)


sql <- "SELECT * FROM fruits;SELECT * FROM test"
res1 <- dbSendQuery(con,sql)
dbFetch(res1, n = -1)
if (dbMoreResults(con)) {
  res2 <- dbNextResult(con)
  dbFetch(res2, n = -1)
}
dbListResults(con)
dbClearResult(res1)
dbClearResult(res2)
dbDisconnect(con1)


con <-dbConnect(dbDriver("MySQL"), dbname = "rmysql", user="root",password="123", host="127.0.0.1", port=3306,client.flag= CLIENT_MULTI_STATEMENTS)
dbListTables(con)
dbReadTable(con,"t_user")
dbWriteTable(con,"fruits",fruits,,append=T,encoding="utf8")
dbReadTable(con,"fruits")

dbGetQuery(con,"SELECT * FROM fruits limit 8")


summary(MySQL(), verbose = TRUE)





