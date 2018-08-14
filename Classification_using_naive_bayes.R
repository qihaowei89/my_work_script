# 
# library(RWeka) # install.packages("RWeka")
# 
# library(rJava) # install.packages("rJava")


install.packages("slam")
library(tm) # install.packages("tm")
library(readr)
library(tidyr) # install.packages("tidyr")
library(magrittr)
library(SnowballC) # install.packages("SnowballC")
library(wordcloud) # install.packages("wordcloud") 
sms= "/home/qihao/workdir/machine_learning/SMSSpamCollection.txt"
sms_raw = separate(read_table(sms,col_names = F),col = X1,into =c("type","text"),sep = "\t")

# sms_raw$type %>% summary()

head(sms_raw)
sms_raw$type <- factor(sms_raw$type)
table(sms_raw$type)

sms_corpus <- VCorpus(VectorSource(sms_raw$text))

inspect(sms_corpus[1:2])

as.character(sms_corpus[[1]])
lapply(sms_corpus[1:2], as.character)


sms_corpus_clean <- tm_map(sms_corpus,content_transformer(tolower))

# as.character(sms_corpus[[1]])

# as.character(sms_corpus_clean[[1]])
sms_corpus_clean <- tm_map(sms_corpus_clean, removeNumbers)
# getTransformations()
# ?stopwords()
sms_corpus_clean <- tm_map(sms_corpus_clean, removeWords, stopwords())
sms_corpus_clean <- tm_map(sms_corpus_clean, removePunctuation)
replacePunctuation <- function(x) { gsub("[[:punct:]]+", " ", x)} # [[:punct:]] 任何标点符号


# wordStem(c("learn", "learned", "learning", "learns"))
sms_corpus_clean <- tm_map(sms_corpus_clean, stemDocument)
sms_corpus_clean <- tm_map(sms_corpus_clean, stripWhitespace)

sms_dtm  <- DocumentTermMatrix(sms_corpus_clean)
sms_dtm2 <- DocumentTermMatrix(sms_corpus, control = list(
  tolower = TRUE,
  removeNumbers = TRUE,
  stopwords = TRUE,
  removePunctuation = TRUE,
  stemming = TRUE
))

sms_dtm3 <- DocumentTermMatrix(sms_corpus, control = list(
  tolower = TRUE,
  removeNumbers = TRUE,
  stopwords = function(x) { removeWords(x, stopwords()) },
  removePunctuation = TRUE,
  stemming = TRUE
))
sms_dtm
sms_dtm2
sms_dtm3

sms_dtm_train <- sms_dtm[1:4180, ]
sms_dtm_test  <- sms_dtm[4181:5574, ]

sms_train_labels <- sms_raw[1:4180, ]$type
sms_test_labels  <- sms_raw[4181:5574, ]$type


prop.table(table(sms_train_labels))
prop.table(table(sms_test_labels))

library(Cairo)


wordcloud(sms_corpus_clean,min.freq = 50, random.order = FALSE)

spam <- subset(sms_raw, type == "spam")
ham <- subset(sms_raw, type == "ham")
table(sms_raw$type)


wordcloud(spam$text, max.words = 40, scale = c(3, 0.5))
wordcloud(ham$text, max.words = 40, scale = c(3, 0.5))





#Chapter 5 Example – identifying risky bank loans using C5.0 decision trees
library(readr)
#Step1 - collecting data 
#Step2 - exploring and preparing the data
credit <- read.csv("/home/wqh/workdir/machine_learning_data/credit.csv")

credit$default[credit$default ==1] <- "No"
credit$default[credit$default ==2] <- "Yes"
credit$default <- as.factor(credit$default)
head(credit)
str(credit)
credit$checking_balance %>% table()
credit[which(credit_train$checking_balance == "> 200 DM"),]$default %>% table
credit[which(credit_train$checking_balance == "unknown"),]$default %>% table
credit[which(credit_train$checking_balance == "1 - 200 DM"),]$default %>% table
credit[which(credit_train$checking_balance == "< 0 DM"),]$default %>% table

table(credit$checking_balance)
table(credit$savings_balance)
table(credit$default)

#Data preparation – creating random training and test datasets
set.seed(123)
train_sample <- sample(1000, 900)
str(train_sample)
credit_train = credit[train_sample,]
credit_test = credit[-train_sample,]

prop.table(table(credit_train$default))
prop.table(table(credit_test$default))

#Step 3 – training a model on the data
# install.packages("C50")
library(C50)
?C5.0Control
credit_train %>% head()

credit_model <- C5.0(credit_train[-17], credit_train$default)

summary(credit_model )

credit_pred <- predict(credit_model, credit_test)
library(gmodels) #install.packages("gmodels")
CrossTable(credit_test$default, credit_pred, prop.chisq = FALSE, prop.c = FALSE, prop.r = FALSE, dnn = c('actual default', 'predicted default'))



#Step 5 – improving model performance
credit_boost10 <- C5.0(credit_train[-17], credit_train$default, trials = 10)
summary(credit_boost10 )

#Boosting the accuracy of decision trees
credit_boost_pred10 <- predict(credit_boost10, credit_test)
CrossTable(credit_test$default, credit_boost_pred10, prop.chisq = FALSE, prop.c = FALSE, prop.r = FALSE, dnn = c('actual default', 'predicted default'))

#Making mistakes more costlier than others
matrix_dimensions <- list(c("No", "Yes"), c("No", "Yes"))
names(matrix_dimensions) <- c("predicted", "actual")
matrix_dimensions


error_cost <- matrix(c(0, 1, 4, 0), nrow = 2, dimnames = matrix_dimensions)

credit_cost <- C5.0(credit_train[-17], credit_train$default,costs = error_cost)

credit_cost_pred <- predict(credit_cost, credit_test)
CrossTable(credit_test$default, credit_cost_pred, prop.chisq = FALSE, prop.c = FALSE, prop.r = FALSE, dnn = c('actual default', 'predicted default'))




