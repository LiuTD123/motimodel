rm(list = ls())

dataset = "train"
datatype = "DFS"
foldpath <- paste0("E:/workdir/08multimodel/09riskbaseline/",datatype)

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

library(tableone)
library(readxl)
# library(xlsx)
# BiocManager::install("xlsx")

# train <- read_excel("E:/workdir/08multimodel/data_train.xlsx")
# test <- read_excel("E:/workdir/08multimodel/data_test.xlsx")
# external <- read_excel("E:/workdir/08multimodel/data_external.xlsx")
# test$Group <- "Test"
# train$Group <- "Train"
# external$Group <- "External"

load(paste0("../../04auccompare/COX/",datatype,"/trainROC_all.RData"))
roc_train <- ROC_all
roc_train$Group <- "Train"
load(paste0("../../04auccompare/COX/",datatype,"/testROC_all.RData"))
roc_test <- ROC_all
roc_test$Group <- "Test"
load(paste0("../../04auccompare/COX/",datatype,"/externalROC_all.RData"))
roc_external <- ROC_all
roc_external$Group <- "External"

data <- rbind(roc_train,roc_test,roc_external)

colnames(data)
## 指定需要分析的变量
all_variable <- colnames(data)[c(1:3,5:14)]

## 指定分类变量，否则分类数据也将以定量数据进行展示
factor_variable <- colnames(data)[c(2,5,7:11)]

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
write.csv(table_statistics, '01.table_statistics_dfs.csv')
save.image("baseline.Rdata")
