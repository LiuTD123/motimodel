rm(list = ls())
setwd("E:/workdir/08multimodel/04auccompare/")

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

# 导入风险评分
load("../results/OS/riskscore_superpc.RData")
coxsig <- read.csv("../03cox/cox_sigfeature.csv",header = F)
coxsig <- coxsig$V1
# ----------------训练集--------------------
data_train <- read_excel("../data_train.xlsx")

# data_train <- data_train[,-c(1,2,5,6:9,30:43)]
# colnames(data_train)[c(6,8,9)] <- c("riskscore","OS.time","OS")

ROC_all <- data_train[,c("id","OS_risk","censor","OS")]
colnames(ROC_all) <- c("ID","DL_RS","fustat","futime")
trainriskscore$ID <- as.numeric(trainriskscore$ID)
trainriskscore <- trainriskscore[,c("ID","RS")]

ROC_all <- merge(trainriskscore,ROC_all)

# --------------cp--------------------
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(coxsig, collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_train)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all, newdata = data_train, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$cp_rs <- risk_scores

# -----------ajcc-------------
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(c("T_stage","N_stage","stage"), collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_train)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all, newdata = data_train, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$ajcc_rs <- risk_scores

save(ROC_all,file = "train_risk.RData")
load("train_risk.RData")
# --------------------------------------------------------------------------------
# ----------timeROC----------------
dataspan = seq(360,180*6,180)

library(timeROC)
ROC_merge <- timeROC(T=ROC_all$futime,   
                   delta=ROC_all$fustat,   
                   marker=ROC_all$RS,   
                   cause=1,                #阳性结局指标数值
                   weighting="marginal",   #计算方法，默认为marginal
                   times=dataspan,       #时间点，选取1年，3年和5年的生存率
                   iid=TRUE)
ROC_DL <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$DL_RS,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_cp <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$cp_rs,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_ajcc <- timeROC(T=ROC_all$futime,   
                    delta=ROC_all$fustat,   
                    marker=ROC_all$ajcc_rs,   
                    cause=1,                #阳性结局指标数值
                    weighting="marginal",   #计算方法，默认为marginal
                    times=dataspan,       #时间点，选取1年，3年和5年的生存率
                    iid=TRUE)
# ci = confint(ROC_DL, level = 0.95)$CI_AUC

# --------时间推移----------
pdf("train.pdf",w=7,h=5)
par(mai=c(1,1,1,2.5))
plotAUCcurve(ROC_merge, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
plotAUCcurve(ROC_cp, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#4C8045",add=TRUE)
plotAUCcurve(ROC_ajcc, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#5D90BA",add=TRUE)
# legend("bottomright",c("Multi","DL","CP","AJCC"),col=c("#EB4B17","#2775AB","#4C8045","#5D90BA"),lty=1,lwd=2)
legend(x= 1150,y= 0.9,xpd = TRUE,
       c(paste0("Merged AUC=",round(mean(ROC_merge[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("DL AUC=",round(mean(ROC_DL[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("CP AUC=",round(mean(ROC_cp[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("AJCC AUC=",round(mean(ROC_ajcc[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,2]),2),")\n")),
       col=c("#EB4B17", "#2775AB", '#4C8045',"#5D90BA"),
       lty=1, lwd=2, cex = 0.7)
dev.off()
compare(ROC_DL,ROC_all,adjusted=TRUE)

# mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
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

# ------------------------验证集-----------------------
# ----------------训练集--------------------
data_test <- read_excel("../data_test.xlsx")

ROC_all <- data_test[,c("id","OS_risk","censor","OS")]
colnames(ROC_all) <- c("ID","DL_RS","fustat","futime")
testriskscore$ID <- as.numeric(testriskscore$ID)
testriskscore <- testriskscore[,c("ID","RS")]

ROC_all <- merge(testriskscore,ROC_all)

# --------------cp--------------------
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(coxsig, collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_test)

risk_scores <- predict(cox_model_all, newdata = data_test, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$cp_rs <- risk_scores

# -----------ajcc-------------
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(c("T_stage","N_stage","stage"), collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_test)

risk_scores <- predict(cox_model_all, newdata = data_test, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$ajcc_rs <- risk_scores

save(ROC_all,file = "test_risk.RData")
# --------------------------------------------------------------------------------
ROC_merge <- timeROC(T=ROC_all$futime,   
                     delta=ROC_all$fustat,   
                     marker=ROC_all$RS,   
                     cause=1,                #阳性结局指标数值
                     weighting="marginal",   #计算方法，默认为marginal
                     times=dataspan,       #时间点，选取1年，3年和5年的生存率
                     iid=TRUE)
ROC_DL <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$DL_RS,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_cp <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$cp_rs,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_ajcc <- timeROC(T=ROC_all$futime,   
                    delta=ROC_all$fustat,   
                    marker=ROC_all$ajcc_rs,   
                    cause=1,                #阳性结局指标数值
                    weighting="marginal",   #计算方法，默认为marginal
                    times=dataspan,       #时间点，选取1年，3年和5年的生存率
                    iid=TRUE)
# ci = confint(ROC_DL, level = 0.95)$CI_AUC

# --------时间推移----------
pdf("test.pdf",w=7,h=5)
par(mai=c(1,1,1,2.5))
plotAUCcurve(ROC_merge, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
plotAUCcurve(ROC_cp, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#4C8045",add=TRUE)
plotAUCcurve(ROC_ajcc, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#5D90BA",add=TRUE)
# legend("bottomright",c("Multi","DL","CP","AJCC"),col=c("#EB4B17","#2775AB","#4C8045","#5D90BA"),lty=1,lwd=2)
legend(x= 1150,y= 0.9,xpd = TRUE,
       c(paste0("Merged AUC=",round(mean(ROC_merge[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("DL AUC=",round(mean(ROC_DL[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("CP AUC=",round(mean(ROC_cp[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("AJCC AUC=",round(mean(ROC_ajcc[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,2]),2),")\n")),
       col=c("#EB4B17", "#2775AB", '#4C8045',"#5D90BA"),
       lty=1, lwd=2, cex = 0.7)
dev.off()
compare(ROC_DL,ROC_all,adjusted=TRUE)

# mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
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

# ------------------------验证集-----------------------
# ----------------外部--------------------
data_test <- read_excel("../data_all_external.xlsx")
ROC_all <- data_test[,c("id","OS_risk","censor","OS")]
colnames(ROC_all) <- c("ID","DL_RS","fustat","futime")
# externalriskscore$ID <- as.numeric(externalriskscore$ID)
externalriskscore <- externalriskscore[,c("ID","RS")]

ROC_all <- merge(externalriskscore,ROC_all)

# --------------cp--------------------
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(coxsig, collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_test)

risk_scores <- predict(cox_model_all, newdata = data_test, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$cp_rs <- risk_scores

# -----------ajcc-------------
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(c("T_stage","N_stage","stage"), collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_test)

risk_scores <- predict(cox_model_all, newdata = data_test, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$ajcc_rs <- risk_scores

save(ROC_all,file = "external_risk.RData")
# --------------------------------------------------------------------------------
ROC_merge <- timeROC(T=ROC_all$futime,   
                     delta=ROC_all$fustat,   
                     marker=ROC_all$RS,   
                     cause=1,                #阳性结局指标数值
                     weighting="marginal",   #计算方法，默认为marginal
                     times=dataspan,       #时间点，选取1年，3年和5年的生存率
                     iid=TRUE)
ROC_DL <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$DL_RS,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_cp <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$cp_rs,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_ajcc <- timeROC(T=ROC_all$futime,   
                    delta=ROC_all$fustat,   
                    marker=ROC_all$ajcc_rs,   
                    cause=1,                #阳性结局指标数值
                    weighting="marginal",   #计算方法，默认为marginal
                    times=dataspan,       #时间点，选取1年，3年和5年的生存率
                    iid=TRUE)
# ci = confint(ROC_DL, level = 0.95)$CI_AUC

# --------时间推移----------
pdf("external.pdf",w=7,h=5)
par(mai=c(1,1,1,2.5))
plotAUCcurve(ROC_merge, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
plotAUCcurve(ROC_cp, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#4C8045",add=TRUE)
plotAUCcurve(ROC_ajcc, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#5D90BA",add=TRUE)
# legend("bottomright",c("Multi","DL","CP","AJCC"),col=c("#EB4B17","#2775AB","#4C8045","#5D90BA"),lty=1,lwd=2)
legend(x= 1150,y= 0.9,xpd = TRUE,
       c(paste0("Merged AUC=",round(mean(ROC_merge[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("DL AUC=",round(mean(ROC_DL[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("CP AUC=",round(mean(ROC_cp[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("AJCC AUC=",round(mean(ROC_ajcc[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,2]),2),")\n")),
       col=c("#EB4B17", "#2775AB", '#4C8045',"#5D90BA"),
       lty=1, lwd=2, cex = 0.7)
dev.off()
compare(ROC_DL,ROC_all,adjusted=TRUE)

compare(ROC_DL_test,ROC_all_test)

# ----------------DFS--------------------
# ----------------训练集--------------------
# 导入风险评分
load("../results/DFS/riskscore_superpc.RData")
coxsig <- read.csv("../03cox/cox_sigfeature.csv",header = F)
coxsig <- coxsig$V1
# ----------------训练集--------------------
data_train <- read_excel("../data_train.xlsx")

# data_train <- data_train[,-c(1,2,5,6:9,30:43)]
# colnames(data_train)[c(6,8,9)] <- c("riskscore","OS.time","OS")

ROC_all <- data_train[,c("id","DFS_risk","recurrence","DFS")]
colnames(ROC_all) <- c("ID","DL_RS","fustat","futime")
trainriskscore$ID <- as.numeric(trainriskscore$ID)
trainriskscore <- trainriskscore[,c("ID","RS")]

ROC_all <- merge(trainriskscore,ROC_all)

# --------------cp--------------------
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(coxsig, collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_train)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all, newdata = data_train, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$cp_rs <- risk_scores

# -----------ajcc-------------
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(c("T_stage","N_stage","stage"), collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_train)
# 
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all, newdata = data_train, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$ajcc_rs <- risk_scores

save(ROC_all,file = "train_risk.RData")
load("train_risk.RData")
# --------------------------------------------------------------------------------
# ----------timeROC----------------
dataspan = seq(360,180*6,180)

library(timeROC)
ROC_merge <- timeROC(T=ROC_all$futime,   
                     delta=ROC_all$fustat,   
                     marker=ROC_all$RS,   
                     cause=1,                #阳性结局指标数值
                     weighting="marginal",   #计算方法，默认为marginal
                     times=dataspan,       #时间点，选取1年，3年和5年的生存率
                     iid=TRUE)
ROC_DL <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$DL_RS,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_cp <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$cp_rs,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_ajcc <- timeROC(T=ROC_all$futime,   
                    delta=ROC_all$fustat,   
                    marker=ROC_all$ajcc_rs,   
                    cause=1,                #阳性结局指标数值
                    weighting="marginal",   #计算方法，默认为marginal
                    times=dataspan,       #时间点，选取1年，3年和5年的生存率
                    iid=TRUE)
# ci = confint(ROC_DL, level = 0.95)$CI_AUC

# --------时间推移----------
pdf("train.pdf",w=7,h=5)
par(mai=c(1,1,1,2.5))
plotAUCcurve(ROC_merge, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
plotAUCcurve(ROC_cp, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#4C8045",add=TRUE)
plotAUCcurve(ROC_ajcc, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#5D90BA",add=TRUE)
# legend("bottomright",c("Multi","DL","CP","AJCC"),col=c("#EB4B17","#2775AB","#4C8045","#5D90BA"),lty=1,lwd=2)
legend(x= 1150,y= 0.9,xpd = TRUE,
       c(paste0("Merged AUC=",round(mean(ROC_merge[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("DL AUC=",round(mean(ROC_DL[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("CP AUC=",round(mean(ROC_cp[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("AJCC AUC=",round(mean(ROC_ajcc[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,2]),2),")\n")),
       col=c("#EB4B17", "#2775AB", '#4C8045',"#5D90BA"),
       lty=1, lwd=2, cex = 0.7)
dev.off()
compare(ROC_DL,ROC_all,adjusted=TRUE)

# mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
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

# ------------------------验证集-----------------------
# ----------------训练集--------------------
data_test <- read_excel("../data_test.xlsx")

ROC_all <- data_test[,c("id","DFS_risk","recurrence","DFS")]
colnames(ROC_all) <- c("ID","DL_RS","fustat","futime")
testriskscore$ID <- as.numeric(testriskscore$ID)
testriskscore <- testriskscore[,c("ID","RS")]

ROC_all <- merge(testriskscore,ROC_all)

# --------------cp--------------------
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(coxsig, collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_test)

risk_scores <- predict(cox_model_all, newdata = data_test, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$cp_rs <- risk_scores

# -----------ajcc-------------
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(c("T_stage","N_stage","stage"), collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_test)

risk_scores <- predict(cox_model_all, newdata = data_test, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$ajcc_rs <- risk_scores

save(ROC_all,file = "test_risk.RData")
# --------------------------------------------------------------------------------
ROC_merge <- timeROC(T=ROC_all$futime,   
                     delta=ROC_all$fustat,   
                     marker=ROC_all$RS,   
                     cause=1,                #阳性结局指标数值
                     weighting="marginal",   #计算方法，默认为marginal
                     times=dataspan,       #时间点，选取1年，3年和5年的生存率
                     iid=TRUE)
ROC_DL <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$DL_RS,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_cp <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$cp_rs,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_ajcc <- timeROC(T=ROC_all$futime,   
                    delta=ROC_all$fustat,   
                    marker=ROC_all$ajcc_rs,   
                    cause=1,                #阳性结局指标数值
                    weighting="marginal",   #计算方法，默认为marginal
                    times=dataspan,       #时间点，选取1年，3年和5年的生存率
                    iid=TRUE)
# ci = confint(ROC_DL, level = 0.95)$CI_AUC

# --------时间推移----------
pdf("test.pdf",w=7,h=5)
par(mai=c(1,1,1,2.5))
plotAUCcurve(ROC_merge, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
plotAUCcurve(ROC_cp, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#4C8045",add=TRUE)
plotAUCcurve(ROC_ajcc, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#5D90BA",add=TRUE)
# legend("bottomright",c("Multi","DL","CP","AJCC"),col=c("#EB4B17","#2775AB","#4C8045","#5D90BA"),lty=1,lwd=2)
legend(x= 1150,y= 0.9,xpd = TRUE,
       c(paste0("Merged AUC=",round(mean(ROC_merge[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("DL AUC=",round(mean(ROC_DL[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("CP AUC=",round(mean(ROC_cp[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("AJCC AUC=",round(mean(ROC_ajcc[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,2]),2),")\n")),
       col=c("#EB4B17", "#2775AB", '#4C8045',"#5D90BA"),
       lty=1, lwd=2, cex = 0.7)
dev.off()
compare(ROC_DL,ROC_all,adjusted=TRUE)

# mycol=rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
# 
# title = paste0("Cli and Gene Risk\n Year = ",cutoff_1,",",cutoff_2,",",cutoff_3)
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

# ------------------------验证集-----------------------
# ----------------外部--------------------
data_test <- read_excel("../data_all_external.xlsx")
ROC_all <- data_test[,c("id","DFS_risk","recurrence","DFS")]
colnames(ROC_all) <- c("ID","DL_RS","fustat","futime")
# externalriskscore$ID <- as.numeric(externalriskscore$ID)
externalriskscore <- externalriskscore[,c("ID","RS")]

ROC_all <- merge(externalriskscore,ROC_all)

# --------------cp--------------------
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(coxsig, collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_test)

risk_scores <- predict(cox_model_all, newdata = data_test, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$cp_rs <- risk_scores

# -----------ajcc-------------
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(c("T_stage","N_stage","stage"), collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = data_test)

risk_scores <- predict(cox_model_all, newdata = data_test, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$ajcc_rs <- risk_scores

save(ROC_all,file = "external_risk.RData")
# --------------------------------------------------------------------------------
ROC_merge <- timeROC(T=ROC_all$futime,   
                     delta=ROC_all$fustat,   
                     marker=ROC_all$RS,   
                     cause=1,                #阳性结局指标数值
                     weighting="marginal",   #计算方法，默认为marginal
                     times=dataspan,       #时间点，选取1年，3年和5年的生存率
                     iid=TRUE)
ROC_DL <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$DL_RS,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_cp <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$cp_rs,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_ajcc <- timeROC(T=ROC_all$futime,   
                    delta=ROC_all$fustat,   
                    marker=ROC_all$ajcc_rs,   
                    cause=1,                #阳性结局指标数值
                    weighting="marginal",   #计算方法，默认为marginal
                    times=dataspan,       #时间点，选取1年，3年和5年的生存率
                    iid=TRUE)
# ci = confint(ROC_DL, level = 0.95)$CI_AUC

# --------时间推移----------
pdf("external.pdf",w=7,h=5)
par(mai=c(1,1,1,2.5))
plotAUCcurve(ROC_merge, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
plotAUCcurve(ROC_cp, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#4C8045",add=TRUE)
plotAUCcurve(ROC_ajcc, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#5D90BA",add=TRUE)
# legend("bottomright",c("Multi","DL","CP","AJCC"),col=c("#EB4B17","#2775AB","#4C8045","#5D90BA"),lty=1,lwd=2)
legend(x= 1150,y= 0.9,xpd = TRUE,
       c(paste0("Merged AUC=",round(mean(ROC_merge[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("DL AUC=",round(mean(ROC_DL[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("CP AUC=",round(mean(ROC_cp[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,2]),2),")\n"),
         paste0("AJCC AUC=",round(mean(ROC_ajcc[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,1]),2),"-",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,2]),2),")\n")),
       col=c("#EB4B17", "#2775AB", '#4C8045',"#5D90BA"),
       lty=1, lwd=2, cex = 0.7)
dev.off()
compare(ROC_DL,ROC_all,adjusted=TRUE)

compare(ROC_DL_test,ROC_all_test)