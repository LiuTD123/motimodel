rm(list = ls())
setwd("E:/workdir/08multimodel/02multi/")

library(gghalves)
library(tidyverse)
library(org.Hs.eg.db)
library(GOSemSim)
library(reshape2)
library(corrplot)
library(ggcorrplot)
library(corrplot)
library(graphics)
library(WGCNA)
library(rstatix)
library(IOBR)
library(reshape2)
library(stats)
library(dplyr)
library(psych)
library(Hmisc)
library(tidyr)
library(survivalROC)
library(survminer)
library(regplot)
library(survival)
library(rms)
library(forestplot)
library(tidyr)
library(survival)
library(patchwork)
library(magrittr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

library(readxl)
# ----------------训练集--------------------
data_train <- read_excel("../data_train2.xlsx")

# data_train <- data_train[,-c(1,2,5,6:9,30:43)]
colnames(data_train)[c(6,8,9,10)] <- c("riskscore","OS.time","OS","ID")

# data_train <- data_train[,c(3,2,1,4:22)]
features <- read.csv("../results/modelgenes.csv")

data_train <- data_train[,c("riskscore","OS","OS.time",features$x)]

rt_all <- data_train

ROC_DL <- data_train[,c("riskscore","OS","OS.time")]

surv_obj_all <- Surv(time = rt_all$OS.time, event = rt_all$OS)
# 
# # 拟合Cox模型
cox_model_all <- coxph(surv_obj_all ~riskscore+ pathologicalgrade+
                         LVI+
                         Nstage+
                         CA199+
                         tumorsize, 
                       data = rt_all)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all, newdata = rt_all, type = "risk")

ROC_all <- rt_all
# # 将风险评分添加到原始数据框中
ROC_all$riskScore <- risk_scores

# --------------------------------------------------------------------------------
ROC_all = subset(ROC_all, select = c(OS, OS.time, riskScore))
colnames(ROC_all) <- c("fustat","futime","riskscore")
colnames(ROC_DL) <- c("riskscore","fustat","futime")

# ----------timeROC----------------
dataspan = seq(180,180*10,180)

library(timeROC)
ROC_all <- timeROC(T=ROC_all$futime,   
                delta=ROC_all$fustat,   
                marker=ROC_all$riskscore,   
                cause=1,                #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=dataspan,       #时间点，选取1年，3年和5年的生存率
                iid=TRUE)
# cutoff_1 <- 1*365
# cutoff_2 <- 3*365
# cutoff_3 <- 5*365
ROC_DL <- timeROC(T=ROC_DL$futime,   
                   delta=ROC_DL$fustat,   
                   marker=ROC_DL$riskscore,   
                   cause=1,                #阳性结局指标数值
                   weighting="marginal",   #计算方法，默认为marginal
                   times=dataspan,       #时间点，选取1年，3年和5年的生存率
                   # times=c(cutoff_1,cutoff_2,cutoff_3),
                   iid=TRUE)
ci = confint(ROC_DL, level = 0.95)$CI_AUC
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
# 
# pdf(file = paste0("01.ModelROC_train_specific.pdf"),width = 5,height = 5)
plot(ROC_DL,
     time=cutoff_1, col="#EB4B17", lwd=2, title="COAD-train")  #time是时间点，col是线条颜色
plot(ROC_DL,
     time=cutoff_2, col="#2775AB", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC_DL,
     time=cutoff_3, col='#4C8045', add=TRUE, lwd=2)
# title(title)
#添加标签信息
legend("bottomright",
       c(paste0("AUC at ",cutoff_1," year: ",round(ROC_DL[["AUC"]][1],3),"\nCI(",ci[1,][1],"-",ci[1,][2],")\n"),
         paste0("AUC at ",cutoff_2," year: ",round(ROC_DL[["AUC"]][2],3),"\nCI(",ci[2,][1],"-",ci[2,][2],")\n"),
         paste0("AUC at ",cutoff_3," year: ",round(ROC_DL[["AUC"]][3],3),"\nCI(",ci[3,][1],"-",ci[3,][2],")")),
       col=c("#EB4B17", "#2775AB", '#4C8045'),
       lty=1, lwd=2,bty = "n")
# dev.off()


pdf("train.pdf",w=5,h=5)
plotAUCcurve(ROC_all, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
legend("bottomright",c("Multi","DL"),col=c("#EB4B17","#2775AB"),lty=1,lwd=2)
dev.off()
compare(ROC_DL,ROC_all,adjusted=TRUE)

# ------------------------验证集-----------------------
# ----------------训练集--------------------
data_test <- read_excel("../data_test2.xlsx")

# data_train <- data_train[,-c(1,2,5,6:9,30:43)]
colnames(data_test)[c(6,8,9,10)] <- c("riskscore","OS.time","OS","ID")

# data_train <- data_train[,c(3,2,1,4:22)]
features <- read.csv("../results/modelgenes.csv")

data_test <- data_test[,c("riskscore","OS","OS.time",features$x)]

rt_all_test <- data_test

ROC_DL_test <- data_test[,c("riskscore","OS","OS.time")]

surv_obj_all_test <- Surv(time = rt_all_test$OS.time, event = rt_all_test$OS)
# 
# # 拟合Cox模型
cox_model_all_test <- coxph(surv_obj_all_test ~riskscore+ pathologicalgrade+
                         LVI+
                         Nstage+
                         CA199+
                         tumorsize, 
                       data = rt_all_test)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all_test, newdata = rt_all_test, type = "risk")

ROC_all_test <- rt_all_test
# # 将风险评分添加到原始数据框中
ROC_all_test$riskScore <- risk_scores

# --------------------------------------------------------------------------------
ROC_all_test = subset(ROC_all_test, select = c(OS, OS.time, riskScore))
colnames(ROC_all_test)[1:3] <- c("fustat","futime","riskscore")
colnames(ROC_DL_test)[1:3] <- c("riskscore","fustat","futime")

# ----------timeROC----------------
dataspan = seq(0,180*10,180)

library(timeROC)
ROC_all_test <- timeROC(T=ROC_all_test$futime,   
                   delta=ROC_all_test$fustat,   
                   marker=ROC_all_test$riskscore,   
                   cause=1,                #阳性结局指标数值
                   # weighting="aalen",   #计算方法，默认为marginal
                   times=dataspan,       #时间点，选取1年，3年和5年的生存率
                   iid=TRUE
                   # other_markers = as.matrix(data_test[,features$x])
                   )
# cutoff_1 <- 1*365
# cutoff_2 <- 3
# cutoff_3 <- 5
ROC_DL_test <- timeROC(T=ROC_DL_test$futime,   
                  delta=ROC_DL_test$fustat,   
                  marker=ROC_DL_test$riskscore,   
                  cause=1,                #阳性结局指标数值
                  # weighting="aalen",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  # other_markers = as.matrix(data_test[,features$x])
                  iid=TRUE
                  )
# ci = confint(ROC1, level = 0.95)$CI_AUC
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
# 
# pdf(file = paste0("01.ModelROC_train_specific.pdf"),width = 5,height = 5)
plot(ROC_DL_test,
     time=cutoff_1, col="#EB4B17", lwd=2, title="COAD-train")  #time是时间点，col是线条颜色
plot(ROC_DL,
     time=cutoff_2, col="#2775AB", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC_DL,
     time=cutoff_3, col='#4C8045', add=TRUE, lwd=2)
title(title)
#添加标签信息
legend("bottomright",
       c(paste0("AUC at ",cutoff_1," year: ",round(ROC1[["AUC"]][1],3),"\nCI(",ci[1,][1],"-",ci[1,][2],")\n"),
         paste0("AUC at ",cutoff_2," year: ",round(ROC1[["AUC"]][2],3),"\nCI(",ci[2,][1],"-",ci[2,][2],")\n"),
         paste0("AUC at ",cutoff_3," year: ",round(ROC1[["AUC"]][3],3),"\nCI(",ci[3,][1],"-",ci[3,][2],")")),
       col=c("#EB4B17", "#2775AB", '#4C8045'),
       lty=1, lwd=2,bty = "n")
# dev.off()
pdf("internaltest.pdf",w=5,h=5)
plotAUCcurve(ROC_all_test, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL_test, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
legend("bottomright",c("Multi","DL"),col=c("#EB4B17","#2775AB"),lty=1,lwd=2)
dev.off()
compare(ROC_DL_test,ROC_all_test)

# ------------------------验证集-----------------------
# ----------------外部--------------------
data_test <- read_excel("../data_all_external2.xlsx")

# data_train <- data_train[,-c(1,2,5,6:9,30:43)]
colnames(data_test)[c(6,8,9,10)] <- c("riskscore","OS.time","OS","ID")

# data_train <- data_train[,c(3,2,1,4:22)]
features <- read.csv("../results/modelgenes.csv")

data_test <- data_test[,c("riskscore","OS","OS.time",features$x)]

rt_all_test <- data_test

ROC_DL_test <- data_test[,c("riskscore","OS","OS.time")]

surv_obj_all_test <- Surv(time = rt_all_test$OS.time, event = rt_all_test$OS)
# 
# # 拟合Cox模型
cox_model_all_test <- coxph(surv_obj_all_test ~riskscore+ pathologicalgrade+
                              LVI+
                              Nstage+
                              CA199+
                              tumorsize, 
                            data = rt_all_test)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all_test, newdata = rt_all_test, type = "risk")

ROC_all_test <- rt_all_test
# # 将风险评分添加到原始数据框中
ROC_all_test$riskScore <- risk_scores

# --------------------------------------------------------------------------------
ROC_all_test = subset(ROC_all_test, select = c(OS, OS.time, riskScore))
colnames(ROC_all_test)[1:3] <- c("fustat","futime","riskscore")
colnames(ROC_DL_test)[1:3] <- c("riskscore","fustat","futime")

# ----------timeROC----------------
dataspan = seq(0,180*10,180)

library(timeROC)
ROC_all_test <- timeROC(T=ROC_all_test$futime,   
                        delta=ROC_all_test$fustat,   
                        marker=ROC_all_test$riskscore,   
                        cause=1,                #阳性结局指标数值
                        # weighting="aalen",   #计算方法，默认为marginal
                        times=dataspan,       #时间点，选取1年，3年和5年的生存率
                        iid=TRUE
                        # other_markers = as.matrix(data_test[,features$x])
)
# cutoff_1 <- 1*365
# cutoff_2 <- 3
# cutoff_3 <- 5
ROC_DL_test <- timeROC(T=ROC_DL_test$futime,   
                       delta=ROC_DL_test$fustat,   
                       marker=ROC_DL_test$riskscore,   
                       cause=1,                #阳性结局指标数值
                       # weighting="aalen",   #计算方法，默认为marginal
                       times=dataspan,       #时间点，选取1年，3年和5年的生存率
                       # other_markers = as.matrix(data_test[,features$x])
                       iid=TRUE
)
# ci = confint(ROC1, level = 0.95)$CI_AUC
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
# 
# pdf(file = paste0("01.ModelROC_train_specific.pdf"),width = 5,height = 5)
plot(ROC_DL_test,
     time=cutoff_1, col="#EB4B17", lwd=2, title="COAD-train")  #time是时间点，col是线条颜色
plot(ROC_DL,
     time=cutoff_2, col="#2775AB", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC_DL,
     time=cutoff_3, col='#4C8045', add=TRUE, lwd=2)
title(title)
#添加标签信息
legend("bottomright",
       c(paste0("AUC at ",cutoff_1," year: ",round(ROC1[["AUC"]][1],3),"\nCI(",ci[1,][1],"-",ci[1,][2],")\n"),
         paste0("AUC at ",cutoff_2," year: ",round(ROC1[["AUC"]][2],3),"\nCI(",ci[2,][1],"-",ci[2,][2],")\n"),
         paste0("AUC at ",cutoff_3," year: ",round(ROC1[["AUC"]][3],3),"\nCI(",ci[3,][1],"-",ci[3,][2],")")),
       col=c("#EB4B17", "#2775AB", '#4C8045'),
       lty=1, lwd=2,bty = "n")
# dev.off()
pdf("external.pdf",w=5,h=5)
plotAUCcurve(ROC_all_test, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL_test, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
legend("bottomright",c("Multi","DL"),col=c("#EB4B17","#2775AB"),lty=1,lwd=2)
dev.off()

compare(ROC_DL_test,ROC_all_test)

# ----------------DFS--------------------
# ----------------训练集--------------------
data_train <- read_excel("../data_train2.xlsx")

# data_train <- data_train[,-c(1,2,5,6:9,30:43)]
colnames(data_train)[c(2,8,5,10)] <- c("riskscore","OS.time","OS","ID")

# data_train <- data_train[,c(3,2,1,4:22)]
features <- read.csv("../results/modelgenes.csv")

data_train <- data_train[,c("riskscore","OS","OS.time",features$x)]

rt_all <- data_train

ROC_DL <- data_train[,c("riskscore","OS","OS.time")]

surv_obj_all <- Surv(time = rt_all$OS.time, event = rt_all$OS)
# 
# # 拟合Cox模型
cox_model_all <- coxph(surv_obj_all ~riskscore+ pathologicalgrade+
                         LVI+
                         Nstage+
                         CA199+
                         tumorsize, 
                       data = rt_all)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all, newdata = rt_all, type = "risk")

ROC_all <- rt_all
# # 将风险评分添加到原始数据框中
ROC_all$riskScore <- risk_scores

# --------------------------------------------------------------------------------
ROC_all = subset(ROC_all, select = c(OS, OS.time, riskScore))
colnames(ROC_all) <- c("fustat","futime","riskscore")
colnames(ROC_DL) <- c("riskscore","fustat","futime")

# ----------timeROC----------------
dataspan = seq(180,180*10,180)

library(timeROC)
ROC_all <- timeROC(T=ROC_all$futime,   
                   delta=ROC_all$fustat,   
                   marker=ROC_all$riskscore,   
                   cause=1,                #阳性结局指标数值
                   weighting="marginal",   #计算方法，默认为marginal
                   times=dataspan,       #时间点，选取1年，3年和5年的生存率
                   iid=TRUE)
# cutoff_1 <- 1*365
# cutoff_2 <- 3*365
# cutoff_3 <- 5*365
ROC_DL <- timeROC(T=ROC_DL$futime,   
                  delta=ROC_DL$fustat,   
                  marker=ROC_DL$riskscore,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  # times=c(cutoff_1,cutoff_2,cutoff_3),
                  iid=TRUE)
ci = confint(ROC_DL, level = 0.95)$CI_AUC
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
# 
# pdf(file = paste0("01.ModelROC_train_specific.pdf"),width = 5,height = 5)
plot(ROC_DL,
     time=cutoff_1, col="#EB4B17", lwd=2, title="COAD-train")  #time是时间点，col是线条颜色
plot(ROC_DL,
     time=cutoff_2, col="#2775AB", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC_DL,
     time=cutoff_3, col='#4C8045', add=TRUE, lwd=2)
# title(title)
#添加标签信息
legend("bottomright",
       c(paste0("AUC at ",cutoff_1," year: ",round(ROC_DL[["AUC"]][1],3),"\nCI(",ci[1,][1],"-",ci[1,][2],")\n"),
         paste0("AUC at ",cutoff_2," year: ",round(ROC_DL[["AUC"]][2],3),"\nCI(",ci[2,][1],"-",ci[2,][2],")\n"),
         paste0("AUC at ",cutoff_3," year: ",round(ROC_DL[["AUC"]][3],3),"\nCI(",ci[3,][1],"-",ci[3,][2],")")),
       col=c("#EB4B17", "#2775AB", '#4C8045'),
       lty=1, lwd=2,bty = "n")
# dev.off()


pdf("train_DFS.pdf",w=5,h=5)
plotAUCcurve(ROC_all, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
legend("bottomright",c("Multi","DL"),col=c("#EB4B17","#2775AB"),lty=1,lwd=2)
dev.off()
compare(ROC_DL,ROC_all,adjusted=TRUE)

# ------------------------验证集-----------------------
# ----------------训练集--------------------
data_test <- read_excel("../data_test2.xlsx")

# data_train <- data_train[,-c(1,2,5,6:9,30:43)]
colnames(data_test)[c(2,8,5,10)] <- c("riskscore","OS.time","OS","ID")

# data_train <- data_train[,c(3,2,1,4:22)]
features <- read.csv("../results/modelgenes.csv")

data_test <- data_test[,c("riskscore","OS","OS.time",features$x)]

rt_all_test <- data_test

ROC_DL_test <- data_test[,c("riskscore","OS","OS.time")]

surv_obj_all_test <- Surv(time = rt_all_test$OS.time, event = rt_all_test$OS)
# 
# # 拟合Cox模型
cox_model_all_test <- coxph(surv_obj_all_test ~riskscore+ pathologicalgrade+
                              LVI+
                              Nstage+
                              CA199+
                              tumorsize, 
                            data = rt_all_test)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all_test, newdata = rt_all_test, type = "risk")

ROC_all_test <- rt_all_test
# # 将风险评分添加到原始数据框中
ROC_all_test$riskScore <- risk_scores

# --------------------------------------------------------------------------------
ROC_all_test = subset(ROC_all_test, select = c(OS, OS.time, riskScore))
colnames(ROC_all_test)[1:3] <- c("fustat","futime","riskscore")
colnames(ROC_DL_test)[1:3] <- c("riskscore","fustat","futime")

# ----------timeROC----------------
dataspan = seq(0,180*10,180)

library(timeROC)
ROC_all_test <- timeROC(T=ROC_all_test$futime,   
                        delta=ROC_all_test$fustat,   
                        marker=ROC_all_test$riskscore,   
                        cause=1,                #阳性结局指标数值
                        # weighting="aalen",   #计算方法，默认为marginal
                        times=dataspan,       #时间点，选取1年，3年和5年的生存率
                        iid=TRUE
                        # other_markers = as.matrix(data_test[,features$x])
)
# cutoff_1 <- 1*365
# cutoff_2 <- 3
# cutoff_3 <- 5
ROC_DL_test <- timeROC(T=ROC_DL_test$futime,   
                       delta=ROC_DL_test$fustat,   
                       marker=ROC_DL_test$riskscore,   
                       cause=1,                #阳性结局指标数值
                       # weighting="aalen",   #计算方法，默认为marginal
                       times=dataspan,       #时间点，选取1年，3年和5年的生存率
                       # other_markers = as.matrix(data_test[,features$x])
                       iid=TRUE
)
# ci = confint(ROC1, level = 0.95)$CI_AUC
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
# 
# pdf(file = paste0("01.ModelROC_train_specific.pdf"),width = 5,height = 5)
plot(ROC_DL_test,
     time=cutoff_1, col="#EB4B17", lwd=2, title="COAD-train")  #time是时间点，col是线条颜色
plot(ROC_DL,
     time=cutoff_2, col="#2775AB", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC_DL,
     time=cutoff_3, col='#4C8045', add=TRUE, lwd=2)
title(title)
#添加标签信息
legend("bottomright",
       c(paste0("AUC at ",cutoff_1," year: ",round(ROC1[["AUC"]][1],3),"\nCI(",ci[1,][1],"-",ci[1,][2],")\n"),
         paste0("AUC at ",cutoff_2," year: ",round(ROC1[["AUC"]][2],3),"\nCI(",ci[2,][1],"-",ci[2,][2],")\n"),
         paste0("AUC at ",cutoff_3," year: ",round(ROC1[["AUC"]][3],3),"\nCI(",ci[3,][1],"-",ci[3,][2],")")),
       col=c("#EB4B17", "#2775AB", '#4C8045'),
       lty=1, lwd=2,bty = "n")
# dev.off()
pdf("internaltest_DFS.pdf",w=5,h=5)
plotAUCcurve(ROC_all_test, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL_test, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
legend("bottomright",c("Multi","DL"),col=c("#EB4B17","#2775AB"),lty=1,lwd=2)
dev.off()
compare(ROC_DL_test,ROC_all_test)

# ------------------------验证集-----------------------
# ----------------外部--------------------
data_test <- read_excel("../data_all_external2.xlsx")

# data_train <- data_train[,-c(1,2,5,6:9,30:43)]
colnames(data_test)[c(2,8,5,10)] <- c("riskscore","OS.time","OS","ID")

# data_train <- data_train[,c(3,2,1,4:22)]
features <- read.csv("../results/modelgenes.csv")

data_test <- data_test[,c("riskscore","OS","OS.time",features$x)]

rt_all_test <- data_test

ROC_DL_test <- data_test[,c("riskscore","OS","OS.time")]

surv_obj_all_test <- Surv(time = rt_all_test$OS.time, event = rt_all_test$OS)
# 
# # 拟合Cox模型
cox_model_all_test <- coxph(surv_obj_all_test ~riskscore+ pathologicalgrade+
                              LVI+
                              Nstage+
                              CA199+
                              tumorsize, 
                            data = rt_all_test)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all_test, newdata = rt_all_test, type = "risk")

ROC_all_test <- rt_all_test
# # 将风险评分添加到原始数据框中
ROC_all_test$riskScore <- risk_scores

# --------------------------------------------------------------------------------
ROC_all_test = subset(ROC_all_test, select = c(OS, OS.time, riskScore))
colnames(ROC_all_test)[1:3] <- c("fustat","futime","riskscore")
colnames(ROC_DL_test)[1:3] <- c("riskscore","fustat","futime")

# ----------timeROC----------------
dataspan = seq(0,180*10,180)

library(timeROC)
ROC_all_test <- timeROC(T=ROC_all_test$futime,   
                        delta=ROC_all_test$fustat,   
                        marker=ROC_all_test$riskscore,   
                        cause=1,                #阳性结局指标数值
                        # weighting="aalen",   #计算方法，默认为marginal
                        times=dataspan,       #时间点，选取1年，3年和5年的生存率
                        iid=TRUE
                        # other_markers = as.matrix(data_test[,features$x])
)
# cutoff_1 <- 1*365
# cutoff_2 <- 3
# cutoff_3 <- 5
ROC_DL_test <- timeROC(T=ROC_DL_test$futime,   
                       delta=ROC_DL_test$fustat,   
                       marker=ROC_DL_test$riskscore,   
                       cause=1,                #阳性结局指标数值
                       # weighting="aalen",   #计算方法，默认为marginal
                       times=dataspan,       #时间点，选取1年，3年和5年的生存率
                       # other_markers = as.matrix(data_test[,features$x])
                       iid=TRUE
)
# ci = confint(ROC1, level = 0.95)$CI_AUC
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
# 
# pdf(file = paste0("01.ModelROC_train_specific.pdf"),width = 5,height = 5)
plot(ROC_DL_test,
     time=cutoff_1, col="#EB4B17", lwd=2, title="COAD-train")  #time是时间点，col是线条颜色
plot(ROC_DL,
     time=cutoff_2, col="#2775AB", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC_DL,
     time=cutoff_3, col='#4C8045', add=TRUE, lwd=2)
title(title)
#添加标签信息
legend("bottomright",
       c(paste0("AUC at ",cutoff_1," year: ",round(ROC1[["AUC"]][1],3),"\nCI(",ci[1,][1],"-",ci[1,][2],")\n"),
         paste0("AUC at ",cutoff_2," year: ",round(ROC1[["AUC"]][2],3),"\nCI(",ci[2,][1],"-",ci[2,][2],")\n"),
         paste0("AUC at ",cutoff_3," year: ",round(ROC1[["AUC"]][3],3),"\nCI(",ci[3,][1],"-",ci[3,][2],")")),
       col=c("#EB4B17", "#2775AB", '#4C8045'),
       lty=1, lwd=2,bty = "n")
# dev.off()
pdf("external_DFS.pdf",w=5,h=5)
plotAUCcurve(ROC_all_test, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL_test, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
legend("bottomright",c("Multi","DL"),col=c("#EB4B17","#2775AB"),lty=1,lwd=2)
dev.off()

compare(ROC_DL_test,ROC_all_test)
