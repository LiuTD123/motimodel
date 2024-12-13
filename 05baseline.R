# 1. 基线统计
rm(list = ls())
setwd("E:/workdir/08multimodel/05basline/")

library(tableone)
library(readxl)
library(xlsx)
BiocManager::install("xlsx")

train <- read_excel("E:/workdir/08multimodel/data_train.xlsx")
test <- read_excel("E:/workdir/08multimodel/data_test.xlsx")
external <- read_excel("E:/workdir/08multimodel/data_all_external.xlsx")

train <- train[,colnames(train)[1:30]]
train$Group <- "Train"
test <- test[,colnames(train)]
test$Group <- "Test"
external <- external[,colnames(train)[1:30]]
external$Group <- "External"
data <- rbind(train,test,external)

colnames(train)
## 指定需要分析的变量
all_variable <- colnames(train)[c(4,5,8,9,11:17,19:30)]

## 指定分类变量，否则分类数据也将以定量数据进行展示
factor_variable <- colnames(train)[c(5,9,11:17,19,21:26,30)]

# 2. 制作基线表
table_statistics <- CreateTableOne(vars = all_variable, 
                                   strata = 'Group', 
                                   data = data, 
                                   factorVars = factor_variable
)
table_statistics <- print(table_statistics, 
                          showAllLevels = TRUE ## 表示展示分类变量所有分类因子的结果
)
# table_statistics <- as.data.frame(table_statistics)
# write.xlsx(table_statistics, "01.table_statistics.xlsx")
write.csv(table_statistics, '01.table_statistics.csv')
save.image("baseline.Rdata")
