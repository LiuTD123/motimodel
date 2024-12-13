rm(list = ls())

dataset = "train"
datatype = "OS"
foldpath <- paste0("E:/workdir/08multimodel/08cindex/",datatype,"/",dataset)

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)

library(readxl)
library(xtable)
library(tableone)
library(openxlsx)  # 用于读取Excel文件
library(survival)   # 生存分析包
library(survminer)  # 用于绘制生存曲线和计算AUC
library(pROC)       # ROC曲线和AUC计算
library(tableone)   # 生成汇总表
library(ggplot2)
library(boot)
library(MASS)
set.seed(555555)
load(paste0("../../../04auccompare/COX/",datatype,"/",dataset,"ROC_all.RData"))

timespan <- seq(360, 180*6+100, by = 180)

risklist <- c("RS","DL","cp","ajcc")
colors <- c("#EB4B17","#2775AB","#4C8045","#D8D155")

library(compareC)
# 拟合两个Cox模型
RS_cindex <- predict(coxph(Surv(futime, fustat) ~ RS, data = ROC_all))
cp_cindex <- predict(coxph(Surv(futime, fustat) ~ cp, data = ROC_all))

# 比较两个模型
result <- compareC(ROC_all$futime, ROC_all$fustat, RS_cindex, cp_cindex)
pval <- result$pval
if (pval > 0.05){
  p <- "ns"
} else if (0.01 < pval && pval <= 0.05) {
  p <- "*"
} else if (0.001 < pval && pval <= 0.01) {
  p <- "**"
} else {
  p <- "***"
}

pdf(paste0(datatype,"_",dataset,"_cindex.pdf"),w= 7, h =5)
par(mai=c(1,1,1,2.5))
cindexmean <- c()
ciconf <- data.frame()
calculate_ci_boot <- function(data, formula) {
  nboot <- 1000  # 自助抽样的次数
  cindex_boot <- numeric(nboot)
  
  for (b in 1:nboot) {
    # 进行自助抽样
    index <- sample(1:nrow(data), replace = TRUE)
    data_boot <- data[index, ]
    # 计算C - index
    cindex_boot[b] <- survConcordance(formula, data = data_boot)$concordance
  }
  
  # 计算置信区间（这里使用2.5%和97.5%分位数作为95%置信区间）
  ci <- quantile(cindex_boot, c(0.025, 0.975))
  
  return(ci)
}

for (i in 1:length(risklist)){
  riskscore <- risklist[i]
  survform <- as.formula(paste0("Surv(futime, fustat) ~ ",riskscore))
  calculate_cindex <- function(data, time_points) {
    cindex_values <- sapply(time_points, function(time_point) {
      data_subset <- data[data$futime <= time_point, ]
      survConcordance(survform, data = data_subset)$concordance
    })
    return(cindex_values)
  }
  
  cindex <- calculate_cindex(ROC_all, timespan)
  cindexmean <- c(cindexmean,mean(cindex,na.rm=T))
  ci <- calculate_ci_boot(ROC_all, survform)
  ci <- data.frame("lower" = ci[1],"higher"=ci[2])
  ciconf <- rbind(ciconf,ci)
  
  if (riskscore == "cp") {
    y2 = cindex[5]
  }
  
  if (riskscore=="RS"){ 
  y1 <- cindex[5]
  plot(timespan, cindex, lwd=3, type = "l", col = colors[i], 
       ylim = c(0.50, 0.85),
       xlim = c(360,1200),
       xlab = "Time (days)", ylab = "C-index", main = paste0(datatype," C-index in ",dataset), bty = "l", xaxt = "n")
  } else {
    lines(timespan, lwd=3, cindex, col = colors[i])
  }
}
lines(x = c(1100,1100),y = c(y1,y2), lwd=3, col = "black")
text(x = 1120, y = c(y1+y2)/2, labels = paste(p))

legend(x= 1200,y= 0.9,xpd = TRUE,
       c(paste0("Merged C-index=",round(cindexmean[1],3),
                "\n CI(",round(ciconf$lower[1],2),"-",round(ciconf$higher[1],2),")\n"),
         paste0("DL C-index=",round(cindexmean[2],3),
         "\n CI(",round(ciconf$lower[2],2),"-",round(ciconf$higher[2],2),")\n"),
          paste0("CP C-index=",round(cindexmean[3],3),
          "\n CI(",round(ciconf$lower[3],2),"-",round(ciconf$higher[3],2),")\n"),
          paste0("AJCC C-index=",round(cindexmean[4],3),
          "\n CI(",round(ciconf$lower[4],2),"-",round(ciconf$higher[4],2),")\n")),
          col=colors,
                 lty=1, lwd=2, cex = 0.7)
dev.off()
