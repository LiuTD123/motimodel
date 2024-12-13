rm(list = ls())

dataset = "external"
datatype = "DFS"

foldpath <- paste0("E:/workdir/08multimodel/06liex/",dataset,"/",datatype)

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

load(paste0("../../../04auccompare/COX/",datatype,"/",dataset,"ROC_all.RData"))

# -----------------校准曲-线1------------------
# 安装并加载所需的R包
year1 = 1
year2 = 2
year3 = 3
# # 建模并完成计算
set.seed(55555)

varlist <-c("RS","DL","cp","ajcc")

for (var in varlist){
  f1 <- cph(formula =  as.formula(paste("Surv(futime, fustat) ~ ",var)),
            data = ROC_all, x = T,y = T,surv = T, time.inc = 365*year1) # time.inc参数和calibrate的u参数后接天数
  cal1 <- calibrate(f1, cmethod = "KM", method = "boot", u = 365*year1, B = 1000) # m，分组到平均包含 m 个受试者的区间；
  
  f2 <- cph(formula =  as.formula(paste("Surv(futime, fustat) ~ ",var)),
            data = ROC_all, x = T,y = T,surv = T, time.inc = 365*year2) # time.inc参数和calibrate的u参数后接天数
  cal2 <- calibrate(f2, cmethod="KM", method="boot", u = 365*year2, B = 1000)
  
  f3 <- cph(formula =  as.formula(paste("Surv(futime, fustat) ~ ",var)),
            data = ROC_all, x = T,y = T,surv = T, time.inc = 365*year3) # time.inc参数和calibrate的u参数后接天数
  cal3 <- calibrate(f3, cmethod="KM", method="boot", u = 365*year3, B = 1000)
  
  pdf(file = paste0(var,"_nomogram_predicted.pdf"), family = "Times", height = 8, width = 8)
  par(mar=c(5,4,2,3),cex=1.5,family="Times")
  plot(cal1,
       subtitles = F,
       lwd=2,lty=1, ##设置线条形状和尺寸
       errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
       xlab=paste0("Nomogram-Predicted Probability of ",year1,",",year2,",",year3," year ",datatype),#便签
       ylab=paste0("Actual ",year1,",",year2,",",year3," year ",datatype," (proportion)") ,#标签
       col="#00468b",#设置一个颜色
       xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
  plot(cal2,
       add = T,
       subtitles = F,
       lwd=2,lty=1,  ##设置线条宽度和线条类型
       errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
       xlab=paste0("Nomogram-Predicted Probability of ",year1,",",year2,",",year3," year ",datatype),#便签
       ylab=paste0("Actual ",year1,",",year2,",",year3," year ",datatype," (proportion)") ,#标签
       col="#ed0000",#设置一个颜色
       xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
  plot(cal3,
       add = T,
       subtitles = F,
       lwd=2,lty=1, ##设置线条形状和尺寸
       errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
       xlab=paste0("Nomogram-Predicted Probability of ",year1,",",year2,",",year3," year ",datatype),#便签
       ylab=paste0("Actual ",year1,",",year2,",",year3," year ",datatype," (proportion)") ,#标签
       col="#42b540",#设置一个颜色
       xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
  #加上图例
  legend("bottomright", legend=c(paste0(year1,"-year"), paste0(year2,"-year"), paste0(year3,"-year")), 
         col=c("#00468b", "#ed0000", "#42b540"), 
         lwd=2)
  #调整对角线
  abline(0,1,lty=5,lwd=2,col="grey")
  dev.off()
}
