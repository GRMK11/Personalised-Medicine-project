setwd("C:/Users/Krishna/Downloads")
getwd()
train_variants=read.csv("train_variants.csv")
test_variants=read.csv("test_variants.csv")
library("tm")
library("dplyr")
library("tidyr")
library("tibble")
# Importing the text data
train_txt=tibble(text=readLines("Ultimate.csv",skip=1))
train_txt=train_txt %>% separate (text,into = c("ID","clinical_txt"),sep = "\\|\\|",skip=1)
train_txt$ID= gsub("\"", "",train_txt$ID)
test_txt=tibble(text = readLines('Test_text.csv', skip = 1))
test_txt= test_txt %>%  separate (text,into=c("ID","clinical_txt"),sep = "\\|\\|")
test_txt$ID= gsub("\"", "",test_txt$ID)

# for the purpuse of preserving memory in R, we fill time to time delete the data bases which are are not required further
# Text mining : To analyse the text data, we will first make it handy. 
#For this purpuse the prepocessing of text data is required. 
#The preprocessing  should be done in the same way for train text and test text. 
#Hence, the merging of both data set is requied first. 
#For this purpose, we will first merge variant data with text data for both and then will merge train and test.
train=merge(train_variants,train_txt, by = "ID")
test=merge(test_variants,test_txt, by = "ID")

full_data=full_join(train,test)
sum(is.na(full_data$Class))

#Starting preprocessing:converting the text into corpus anddoing feature engineering with tf-idf method.
set.seed(123)
corpus = VCorpus(VectorSource(full_data$clinical_txt))
corpus = tm_map(corpus,content_transformer(stringi::stri_trans_tolower))
corpus = tm_map(corpus,removePunctuation,preserve_intra_word_dashes=T)
corpus = tm_map(corpus, removeWords, stopwords("english"))
corpus = tm_map(corpus, removeNumbers)
# I have created a  my own stopword list
corpus = tm_map(corpus,removeWords,c("author","describe","find","found", "result",
                                     "conclude","analyze","analysis","show","shown","resulted","concluded","described",
                                     "concluded","evaluate","evaluated","discuss","discussed","demonstrate","demonstrated",
                                     "the","this","that","these","those","illustrated","illustrate","list","fig","figure",
                                     "et","al","data","determined","studied","indicated","research","method","determine",
                                     "studies","study","indicate","research","researcher","medical","background","abstract",
                                     "and","but","all","also","are","been","both","can","consider","describe","described",
                                     "declar","determin","did","rt","http\\*"))
corpus = tm_map(corpus, stemDocument, language="english")
corpus = tm_map(corpus, stripWhitespace)

# Converting the corpus into Sparse Matrix for the purpose of analysis
dtm = DocumentTermMatrix(corpus,control = list(weighting = weightTfIdf))
dtm
# No. of terms in the term matrix , which actually represent the no. of columns in the term matrix, which is actually very high. There may exist some very few frequency words.
dtm = removeSparseTerms(dtm, sparse = 0.95) 
dtm

# Identifying most frequest words and their frequence
count_word = colSums(txt_data)
length(count_word)
freq_word=tibble(name = attributes(count_word)$names, count = count_word)
top_n(freq_word,10)
top_n(freq_word,-10)
library("wordcloud")
wordcloud(names(count_word), count_word, min.freq = 100,scale = c(6,.1), colors = brewer.pal(6, 'Dark2'))
dtm
#Dimension R"e1071")
library(e1071)
library("preProcess")

full_data$Gene=as.factor(full_data$Gene)
"class"(full_data$Gene)
full_data$Variation=as.factor(full_data$Variation)
full_data$Variation=as.numeric(full_data$Variation)
full_data$Class=as.factor(full_data$Class)
class(full_data$Gene)
full_data$Gene=as.factor(full_data$Gene)
class(full_data$Gene)
full_data1=subset(full_data,select=c("Class","ID","Gene","Variation"))

newdata=cbind(full_data1,txt_data)
ncol(newdata)
newdata$Class=as.factor(full_data$Class)
newdata$ID=as.factor(newdata$ID)
newdata$ID=as.integer((full_data$ID))
summary(newdata)
str(newdata$Class)
sum(is.na(newdata))

#Splitting the data into train and test
train_newdata = newdata[1:3321,]

nrow(train_newdata)
sum(is.na(train_newdata$Class))
test_newdata = newdata[3322:8989,]
nrow(test_newdata)
test_newdata$Class=NULL
#Splitting train data in  to train and test for building the model
model_train=train_newdata[sample(nrow(train_newdata),2657,replace=F),]
Validation_test=train_newdata[!(1:nrow(train_newdata)) %in% as.numeric(row.names(model_train)),]

#Now we are ready with our data to be train for classification 

#training the model with default parameters
Tmodel=naiveBayes(Class ~., data =model_train)

library("data.table")
# Predicting target variable of validation set through built in model
pred=predict(Tmodel,Validation_test[,2:2674])
#predicting target variable of  validation set with probabilities
pred_prob=predict(Tmodel,Validation_test[,2:2674],type="raw")
pred=as.matrix(pred)
library("MLmetrics")
library("caret")

xtab=ConfusionMatrix(pred,Validation_test$Class)
xtab=as.matrix(xtab)
Accuracy(pred,Validation_test$Class)
ConfusionMatrix(pred,Validation_test$Class)
Accuracy(pred,Validation_test$Class)
RMSE(y_pred=pred_prob,y_true=y_true)
kappa(y_pred=pred_prob,y_true=y_true)
y_true = as.vector(Validation_test$Class)
y_true =c(y_true, c(paste0(seq(1, 9))))
y_true = model.matrix(~ 0 + ., data.frame(as.character(y_true)))
y_true = y_true[-((nrow(y_true)-9+1):nrow(y_true)), ]
Logloss_score=MultiLogLoss(y_pred=pred_prob,y_true=y_true)
Logloss
#Finalmodel to predict our competition testset
Final_model=naiveBayes(Class~.,data=train_newdata)
Final_pred=predict(Final_model,test_newdata[,2:2674])
Final_pred_prob=predict(Final_model,test_newdata[,2:2674],type="raw")
Final_pred=as.numeric(Final_pred)
Final_output_tdata= as.vector(Final_pred)
Final_output_tdata =c(Final_output_tdata, c(paste0(seq(1, 9))))
Final_output_tdata = model.matrix(~ 0 + ., data.frame(as.character(Final_output_tdata)))
Final_output_tdata = Final_output_tdata[-((nrow(Final_output_tdata)-9+1):nrow(Final_output_tdata)), ]
Final_output_tdata=read.csv("Final_output_tdata",header=FALSE)
colnames(Final_output_tdata)[1]= "Class1"
colnames(Final_output_tdata)[2]= "Class2"
colnames(Final_output_tdata)[3]= "Class3"
colnames(Final_output_tdata)[4]= "Class4"
colnames(Final_output_tdata)[5]= "Class5"
colnames(Final_output_tdata)[6]= "Class6"
colnames(Final_output_tdata)[7]= "Class7"
colnames(Final_output_tdata)[8]= "Class8"
colnames(Final_output_tdata)[9]= "Class9"
Final_output_tdata=as.data.frame(Final_output_tdata)
Final_output_tdata=tibble::rowid_to_column(Final_output_tdata, "ID")
#Output file saved as Submissions written back to directory
write.table(data.table(Final_output_tdata), "submission.csv", sep=",", dec=".", quote=FALSE, row.names=FALSE)
