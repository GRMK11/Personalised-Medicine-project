#Adding libraries
library(plyr)
library(dplyr)
library(ggplot2)
library("tidyr")
library("Scale")
#Taking in the Train data
Train = read.csv("../input/training_variants")
setwd("C:/Users/krishna/Downloads")
getwd()
train_variants=read.csv("train_variants.csv")
test_variants=read.csv("test_variants.csv")

test_variants=read.csv("test_variants.csv")
# Importing the text data
train_txt=tibble(text=readLines("Ultimate.csv",skip=1))
train_txt =tibble::text = readLines('Ultimate.csv', skip = 1)
train_txt=train_txt %>% separate (text,into = c("ID","clinical_txt"),sep = "\\|\\|",skip=1)
train_txt$ID= gsub("\"", "",train_txt$ID)
train_txt$ID=as.integer(train_txt$ID)
test_txt=tibble(text = readLines('Test_text.csv', skip = 1))
test_txt= test_txt %>%  separate (text,into=c("ID","clinical_txt"),sep = "\\|\\|")
test_txt$ID= gsub("\"", "",test_txt$ID)
train=merge(train_variants,train_txt, by = "ID")
test=merge(test_variants,test_txt, by = "ID")
train=merge(train_variants,train_txt, by = "ID")
test=merge(test_variants,test_txt, by = "ID")
ggplot(train,aes_string(x=train$Gene,y=train$Class))+
  geom_point(aes_string(colour= train$Class),size=10)+
  theme_bw()+ ylab("Class")+ xlab("Gene") + ggtitle("Main") +
  theme(text=element_text(size=25))+
  scale_y_continuous(breaks=pretty(train$Class,n=10))
#Calculating the frequencies of Variations and sorting in decreasing order
Var= train$Variation
Var.freq = sort(table(Var), decreasing = T)


######################################################################

#Calculating the frequencies of Gene and sorting in decreasing order
Gene = train$Gene
Gene.freq = sort(table(Gene), decreasing = T)

######################################################################

#Calculating the frequencies of class and sorting in decreasing order
class = train$Class
class_freq = data.frame(sort(table(class), decreasing = T))

class_freq[2] = class_freq[1]
class_freq[1] = row.names(class_freq)
class_freq[,2]=as.data.frame(class_freq[,2])
class_freq[,1]=as.data.frame(class_freq[,1])
colnames(class_freq) = c('Class','ClassOccurrances')
class_freq = class_freq[order(class_freq[,1], decreasing = F),]
barplot(class_freq[,2],names.arg = class_freq[,1],
        xlab = 'Class', ylab = 'Frequency', 
        col = c('blue','green','yellow','black','red','maroon','purple','orange','brown'))
library(ggplot2)
ggplot(train, aes(x=Variation, y=Class)) + geom_bar(stat="identity") + 
  labs(x="Gene", y="Class")
######################################################################

#Grouping Gene with Class

Gene_Grp = aggregate(Train$Class, by=list(Train$Gene),FUN= n_distinct)                                                        
Gene_Grp = Gene_Grp[order(Gene_Grp[,2], decreasing = T),]
#View(Gene_Grp)
densityplot(train$Variation, 
            main="Density Plot", 
            xlab="Variation")
#Grouping Variation with Class
library('ggplot2') # visualization
library('ggthemes') # visualization
library('scales') # visualization
library('grid') # visualisation
library('gridExtra') # visualisation
library('corrplot') # visualisation
library('ggfortify') # visualisation
library('ggraph') # visualisation
library('igraph') # visualisation
library('dplyr') # data manipulation
library('readr') # data input
library('tibble') # data wrangling
library('tidyr') # data wrangling
library('stringr') # string manipulation
library('forcats') # factor manipulation
library('tidytext') # text mining
library('SnowballC') # text analysis
library('wordcloud') # test visualisation
Var_Grp = aggregate(train$Class, by=list(train$Variation),FUN= n_distinct)                                                        
Var_Grp = Var_Grp[order(Var_Grp[,2], decreasing = T),]
#View(Var_Grp)
train <- train %>%
  mutate(Gene = factor(Gene),
         Variation = factor(Variation),
         Class = factor(Class))

test <- test %>%
  mutate(Gene = factor(Gene),
         Variation = factor(Variation))

summary(train, maxsum = 9)
glimpse(train)
nrow(train)
## [1] 3321
nrow(test)
## [1] 5668
sum(is.na(train))
## [1] 0
sum(is.na(test))
train %>%
  group_by(Gene) %>%
  summarise(ct = n()) %>%
  arrange(desc(ct))
test %>%
  group_by(Gene) %>%
  summarise(ct = n()) %>%
  arrange(desc(ct))
train %>%
  group_by(Variation) %>%
  summarise(ct = n()) %>%
  arrange(desc(ct))
test %>%
  group_by(Variation) %>%
  summarise(ct = n()) %>%
  arrange(desc(ct))
Individual feature visualisations
#This is the frequency distribution of the most frequent Gene values:
   top_gene <- train %>%
  group_by(Gene) %>%
  summarise(ct = n()) %>%
  filter(ct > 40)
top_gene %>%
  ggplot(aes(reorder(Gene, -ct, FUN = min), ct)) +
  geom_point(size = 4) +
  labs(x = "Gene", y = "Frequency") +
  coord_flip()
top_gene %>%
  ggplot(aes(reorder(Gene, -ct, FUN = min), ct)) +
  geom_point(size = 4) +
  labs(x = "Gene", y = "Frequency") +
  coord_flip()
top_gene_test =test %>%
  group_by(Gene) %>%
  summarise(ct = n()) %>%
  filter(ct > 40)

top_gene_test %>%
  ggplot(aes(reorder(Gene, -ct, FUN = min), ct)) +
  geom_point(size = 4) +
  labs(x = "Gene", y = "Frequency") +
  coord_flip()

#findings
#A relatively small group of Gene levels make up a sizeable part of the feature values in both train and test data.
#The test data has fewer high-frequency Genes.

#These are the most frequent Variations in the train (blue) vs test (red) data; confirming what we already saw by comparing the table data:
######################################################################
train_1=subset(train,select=c("Gene","Variation"))
test_1=subset(test,select=c("Gene","Variation"))
New_1 = train_1 %>% mutate(set = factor("train"))
join_1 = test %>% mutate(set = factor("test"))
New_join = full_join(New_1,join_1)

New_join %>%
  group_by(Variation, set) %>%
  summarise(ct = n()) %>%
  filter(ct > 3) %>%
  ggplot(aes(reorder(Variation, -ct, FUN = median), ct, colour = set)) +
  geom_point(size = 4) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Variation", y = "Frequency")

Here we see how the Class target is distributed in the train data:
  
  train %>%
  ggplot(aes(Class)) +
  geom_bar()
#Feature interactions
#Now we want to examine how the features interact with each other and with the target Class variable.

 #Gene vs Class

#First, we will look at the frequency distribution of the overall most frequent Genes for the different Classes. Note the logarithmic frequency scale.
library("stringr")
train %>%
  filter(Gene %in% str_c(top_gene$Gene)) %>%
  ggplot(aes(Gene)) +
  geom_bar() +
  scale_y_log10() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=7)) +
  facet_wrap(~ Class)
# Classes sorted by Genes (again log counts):
  
  train %>%
  filter(Gene %in% str_c(top_gene$Gene)) %>%
  ggplot(aes(Class)) +
  geom_bar() +
  scale_y_log10() +
  facet_wrap(~ Gene)

