rm(list = ls())

datatype <- "OS"
dataset <- "test"
foldpath <- paste0("E:/workdir/08multimodel/04auccompare/COX/",datatype)

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
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
data_train <- read_excel(paste0("../../../data_",dataset,".xlsx"))

coxsig <- read.csv("../../../03cox/cox_sigfeature.csv",header = F)
coxsig <- coxsig$V1

score <- paste0(datatype,"_risk")
futime <- datatype
if (datatype == "OS"){
  fustat <- "censor"
} else {
  fustat <- "recurrence"
}

ROC_all <- data_train[,c(futime,fustat,score,"id",coxsig,"stage")]
colnames(ROC_all)[1:3] <- c("futime","fustat","DL","ID")

# ---------融合----------
cox_data <- as.formula(paste0('Surv(futime, fustat)~', paste(c("DL",coxsig), collapse = "+")))
# # 拟合Cox模型
cox_model_all <- coxph(cox_data,data = ROC_all)
# # 计算风险评分
# # 使用predict函数和newdata参数来计算新数据集的风险评分
# # 如果您想要计算训练数据集的风险评分，可以使用模型本身作为newdata
risk_scores <- predict(cox_model_all, newdata = ROC_all, type = "risk")

# # 将风险评分添加到原始数据框中
ROC_all$RS <- as.data.frame(risk_scores)$risk_score

# ---------cp----------------
cox_data <- as.formula(paste0('Surv(futime, fustat)~', paste(coxsig, collapse = "+")))
cox_model_all <- coxph(cox_data,data = ROC_all)
risk_scores <- predict(cox_model_all, newdata = ROC_all, type = "risk")
ROC_all$cp <- as.data.frame(risk_scores)$risk_score
# -------- ajcc----------------
cox_data <- as.formula(paste0('Surv(futime, fustat)~', paste(c("T_stage","N_stage","stage"), collapse = "+")))
cox_model_all <- coxph(cox_data,data = ROC_all)
risk_scores <- predict(cox_model_all, newdata = ROC_all, type = "risk")
ROC_all$ajcc <- as.data.frame(risk_scores)$risk_score

save(ROC_all,file = paste0(dataset,"ROC_all.RData"))
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
                  marker=ROC_all$DL,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_cp <- timeROC(T=ROC_all$futime,   
                  delta=ROC_all$fustat,   
                  marker=ROC_all$cp,   
                  cause=1,                #阳性结局指标数值
                  weighting="marginal",   #计算方法，默认为marginal
                  times=dataspan,       #时间点，选取1年，3年和5年的生存率
                  iid=TRUE)
ROC_ajcc <- timeROC(T=ROC_all$futime,   
                    delta=ROC_all$fustat,   
                    marker=ROC_all$ajcc,   
                    cause=1,                #阳性结局指标数值
                    weighting="marginal",   #计算方法，默认为marginal
                    times=dataspan,       #时间点，选取1年，3年和5年的生存率
                    iid=TRUE)
# ci = confint(ROC_DL, level = 0.95)$CI_AUC
pval <- compare(ROC_merge,ROC_cp)
colors <- c("#EB4B17","#2775AB","#4C8045","#D8D155")

pval <- min(pval$p_values_AUC)
if (pval > 0.05){
  p <- "ns"
} else if (0.01 < pval && pval <= 0.05) {
  p <- "*"
} else if (0.001 < pval && pval <= 0.01) {
  p <- "**"
} else {
  p <- "***"
}
y1 = ROC_merge$AUC[5]
y2 = ROC_cp$AUC[5]
# --------时间推移----------
pdf(paste0(dataset,"_",datatype,".pdf"),w=7,h=5)
par(mai=c(1,1,1,2.5))
# plotAUCcurve(ROC_merge, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#EB4B17")
# plotAUCcurve(ROC_DL, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#2775AB",add=TRUE)
# plotAUCcurve(ROC_cp, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#4C8045",add=TRUE)
# plotAUCcurve(ROC_ajcc, FP = 2, conf.int = FALSE, conf.band = FALSE, col = "#D8D155",add=TRUE)

plot(dataspan, ROC_merge$AUC, lwd=3, type = "l", col = "#EB4B17", 
     ylim = c(0.50, 1),
     xlim = c(360,1200),
     xlab = "Time (days)", ylab = "AUC", main = paste0(datatype," AUC in ",dataset), bty = "l", xaxt = "n")
lines(dataspan, lwd=3, ROC_DL$AUC, col = "#2775AB")
lines(dataspan, lwd=3, ROC_cp$AUC, col = "#4C8045")
lines(dataspan, lwd=3, ROC_ajcc$AUC, col = "#D8D155")

legend(x= 1200,y= 0.9,xpd = TRUE,
       c(paste0("Merged AUC=",round(mean(ROC_merge[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,1])/100,2),"-",round(mean(confint(ROC_merge, level = 0.95)$CI_AUC[,2])/100,2),")\n"),
         paste0("DL AUC=",round(mean(ROC_DL[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,1])/100,2),"-",round(mean(confint(ROC_DL, level = 0.95)$CI_AUC[,2]),2)/100,")\n"),
         paste0("CP AUC=",round(mean(ROC_cp[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,1])/100,2),"-",round(mean(confint(ROC_cp, level = 0.95)$CI_AUC[,2]),2)/100,")\n"),
         paste0("AJCC AUC=",round(mean(ROC_ajcc[["AUC"]],na.rm=T),3),
                "\n CI(",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,1])/100,2),"-",round(mean(confint(ROC_ajcc, level = 0.95)$CI_AUC[,2])/100,2),")\n")),
       col=c("#EB4B17", "#2775AB", '#4C8045',"#D8D155"),
       lty=1, lwd=2, cex = 0.7)

lines(x = c(1100,1100),y = c(y1,y2), lwd=3, col = "black")
text(x = 1120, y = c(y1+y2)/2, labels = paste(p))

dev.off()
