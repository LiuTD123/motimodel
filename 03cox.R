rm(list = ls())
setwd("E:/workdir/08multimodel/03cox/")

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
library(writexl)
library(forestplot)

data <- read.xlsx("../data_train.xlsx")
clinical <- c('Age', 'Sex', 'TB', 'Albumin', 'CA199', 'DM', 'HBP', 'BMI', 'location', 'tumor_size', 'imaging_enlarged_LN')
pathological <- c('pathological_grade', 'LVI', 'PNI', 'PFI', 'T_stage', 'N_stage')
CP<- c('Age', 'Sex', 'TB', 'Albumin', 'CA199', 'DM', 'HBP', 'BMI', 'location', 'tumor_size', 'imaging_enlarged_LN',
       'pathological_grade', 'LVI', 'PNI', 'PFI', 'T_stage', 'N_stage')
factor_names <- c('OS_group','DFS_group','SMAD4', 'pathological_grade', 'LVI', 'PNI', 'PFI', 'T_stage', 'N_stage', 'Sex', 'TB', 'Albumin', 'CA199', 'DM', 'HBP', 'location', 'imaging_enlarged_LN')
for(var in factor_names){
  data[[var]] <- factor(data[[var]])
}
# -------------------clinical-----------------
results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

for (var in clinical) {
  formula <- as.formula(paste("Surv(OS, censor) ~", var))
  cox_model <- coxph(formula, data = data)
  summary_cox <- summary(cox_model)
  coef <- summary_cox$coefficients[, "coef"]
  hr <- summary_cox$coefficients[, "exp(coef)"]
  ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[, "Pr(>|z|)"]
  results <- rbind(results, data.frame(Variable = var, Coefficient = coef, HR = hr, 
                                       CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
}

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res2 <- rbind(res) %>% as.data.frame()

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.univariate_cox_os_clinical.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

uni_OS_clinical <- results$Variable[which(results$P_value < 0.05)]
# ---------pathological---------------
results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

for (var in pathological) {
  formula <- as.formula(paste("Surv(OS, censor) ~", var))
  cox_model <- coxph(formula, data = data)
  summary_cox <- summary(cox_model)
  coef <- summary_cox$coefficients[, "coef"]
  hr <- summary_cox$coefficients[, "exp(coef)"]
  ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[, "Pr(>|z|)"]
  if (is.matrix(ci)){
    for (i in 1:nrow(ci)){
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef[i], HR = hr[i], 
                                           CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
    }} else {
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef, HR = hr, 
                                           CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
    }
}

uni_OS_pathological <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1:4,9,5,6,10,7,8),]

res2$Indicator <- c(res2$Indicator[1:4],
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3", 
                    "N Stage\n(N1 Reference)", 
                    "N2", 
                    "N3")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.univariate_cox_os_pathological.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
# ---------cp--------------------------
results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)
for (var in CP) {
  formula <- as.formula(paste("Surv(OS, censor) ~", var))
  cox_model <- coxph(formula, data = data)
  summary_cox <- summary(cox_model)
  coef <- summary_cox$coefficients[, "coef"]
  hr <- summary_cox$coefficients[, "exp(coef)"]
  ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[, "Pr(>|z|)"]
  if (is.matrix(ci)){
    for (i in 1:nrow(ci)){
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef[i], HR = hr[i], 
                                           CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
    }} else {
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef, HR = hr, 
                                           CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
    }
}
uni_OS_CP <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1:15,20,16,17,21,18,19),]

res2$Indicator <- c(res2$Indicator[1:15],
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3", 
                    "N Stage\n(N1 Reference)", 
                    "N2", 
                    "N3")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.univariate_cox_os_CP.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

print(OS_clinical)
print(OS_pathological)
print(OS_CP)

inputclinical <- Reduce(intersect,list(uni_OS_clinical,uni_DFS_clinical,multi_os_clinical,multi_DFS_clinical))
inputpathological <- Reduce(intersect,list(uni_OS_pathological,uni_DFS_pathological,multi_os_pathological,multi_DFS_pathological))

inputclinical <- Reduce(intersect,list(uni_os_clinical,multi_os_clinical))
inputpathological <- Reduce(intersect,list(uni_os_pathological,multi_os_pathological))

input <- c("CA199","tumor_size","pathological_grade","PNI","T_stage","N_stage")
write(input, file = "cox_sigfeature.csv")
# -----单因素画图-------------
# 8,8,13
pdf(file = "07.univariate_cox_dfs_cp.pdf", family = "Times", height = 13, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           # is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,2),TRUE,rep(FALSE,4)),
           # is.summary = c(rep(FALSE,5),TRUE,rep(FALSE,2),TRUE,rep(FALSE,2)),
           is.summary = c(rep(FALSE,16),TRUE,rep(FALSE,2),TRUE,rep(FALSE,2)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           # xticks = c(0, 1, 2, 4, 6, 8), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.4,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily = "Times"),
                          xlab=gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily = "Times")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate",
           clip = c(0,7)) # 垂直于x轴的网格线，对应每个刻度
dev.off()
# ----------------------------
##################DFS--------------
# ------------clinical------------
results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)
for (var in clinical) {
  formula <- as.formula(paste("Surv(DFS, recurrence) ~", var))
  cox_model <- coxph(formula, data = data)
  summary_cox <- summary(cox_model)
  coef <- summary_cox$coefficients[, "coef"]
  hr <- summary_cox$coefficients[, "exp(coef)"]
  ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[, "Pr(>|z|)"]
  if (is.matrix(ci)){
    for (i in 1:nrow(ci)){
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef[i], HR = hr[i], 
                                           CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
    }} else {
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef, HR = hr, 
                                           CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
    }
}
uni_DFS_clinical <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res2 <- rbind(res) %>% as.data.frame()

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.univariate_cox_dfs_clinical.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
# -----------pathological-----------------------
results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)
for (var in pathological) {
  formula <- as.formula(paste("Surv(DFS, recurrence) ~", var))
  cox_model <- coxph(formula, data = data)
  summary_cox <- summary(cox_model)
  coef <- summary_cox$coefficients[, "coef"]
  hr <- summary_cox$coefficients[, "exp(coef)"]
  ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[, "Pr(>|z|)"]
  if (is.matrix(ci)){
    for (i in 1:nrow(ci)){
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef[i], HR = hr[i], 
                                           CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
    }} else {
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef, HR = hr, 
                                           CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
    }
}
uni_DFS_pathological <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1:4,9,5,6,10,7,8),]

res2$Indicator <- c(res2$Indicator[1:4],
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3", 
                    "N Stage\n(N1 Reference)", 
                    "N2", 
                    "N3")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.univariate_cox_dfs_pathological.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
# -------------cp--------------------
results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)
for (var in CP) {
  formula <- as.formula(paste("Surv(DFS, recurrence) ~", var))
  cox_model <- coxph(formula, data = data)
  summary_cox <- summary(cox_model)
  coef <- summary_cox$coefficients[, "coef"]
  hr <- summary_cox$coefficients[, "exp(coef)"]
  ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
  p_value <- summary_cox$coefficients[, "Pr(>|z|)"]
  if (is.matrix(ci)){
    for (i in 1:nrow(ci)){
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef[i], HR = hr[i], 
                                           CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
    }} else {
      results <- rbind(results, data.frame(Variable = var, Coefficient = coef, HR = hr, 
                                           CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
    }
}
uni_DFS_CP <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1:15,20,16,17,21,18,19),]

res2$Indicator <- c(res2$Indicator[1:15],
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3", 
                    "N Stage\n(N1 Reference)", 
                    "N2", 
                    "N3")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.univariate_cox_dfs_CP.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

print(DFS_clinical)
print(DFS_pathological)
print(DFS_CP)

# ------------=============-------------------
##mutli-COX
# ------OS---------------
# -----------clinical----------------
clinical <- c("CA199", "tumor_size")
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(clinical, collapse = "+")))
cox_more <- coxph(cox_data, data = data)

cox_zph <- cox.zph(cox_more)
names(cox_more$coefficients) <- rownames(cox_zph$table)[1:2]

mul_cox_result <- summary(cox_more)$coefficients

summary_cox <- summary(cox_more)

coef <- summary_cox$coefficients[, "coef"]
hr <- summary_cox$coefficients[, "exp(coef)"]
ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
p_value <- summary_cox$coefficients[, "Pr(>|z|)"]

results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

if (is.matrix(ci)){
  for (i in 1:nrow(ci)){
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef[i], HR = hr[i], 
                                         CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
  }} else {
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef, HR = hr, 
                                         CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
}

multi_os_clinical <- results$Variable[which(results$P_value < 0.05)]
  
res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res2 <- rbind(res) %>% as.data.frame()

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.multivariate_cox_os_clinical.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

# ------------------pathological-------------------
pathological <- c("pathological_grade","LVI","PNI","T_stage","N_stage")
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(pathological, collapse = "+")))
cox_more <- coxph(cox_data, data = data)

cox_zph <- cox.zph(cox_more)
names(cox_more$coefficients) <- c(pathological[1:3],
                                  "T2", 
                                  "T3",
                                  "N2", 
                                  "N3")

mul_cox_result <- summary(cox_more)$coefficients

summary_cox <- summary(cox_more)

coef <- summary_cox$coefficients[, "coef"]
hr <- summary_cox$coefficients[, "exp(coef)"]
ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
p_value <- summary_cox$coefficients[, "Pr(>|z|)"]

results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

if (is.matrix(ci)){
  for (i in 1:nrow(ci)){
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef[i], HR = hr[i], 
                                         CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
  }} else {
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef, HR = hr, 
                                         CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
  }

multi_os_pathological <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1:3,9,4,5,8,6,7),]

res2$Indicator <- c(res2$Indicator[1:3],
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3", 
                    "N Stage\n(N1 Reference)", 
                    "N2", 
                    "N3")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.multivariate_cox_os_pathological.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

# --------------cp-----------------
CP <- c("CA199","tumor_size","pathological_grade","LVI","PNI","T_stage","N_stage")
cox_data <- as.formula(paste0('Surv(OS, censor)~', paste(CP, collapse = "+")))
cox_more <- coxph(cox_data, data = data)

cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]
names(cox_more$coefficients) <- c(CP[1:5],
                                  "T2", 
                                  "T3",
                                  "N2", 
                                  "N3")

# cox_formula <- as.formula(paste("Surv(OS, censor)~",
#                                 paste(rownames(cox_table)[cox_table[,3]>0.05],
#                                       collapse = "+")))
# cox_more_2 <- coxph(cox_formula, data = data)
mul_cox_result <- summary(cox_more)$coefficients

summary_cox <- summary(cox_more)

coef <- summary_cox$coefficients[, "coef"]
hr <- summary_cox$coefficients[, "exp(coef)"]
ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
p_value <- summary_cox$coefficients[, "Pr(>|z|)"]

results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

if (is.matrix(ci)){
  for (i in 1:nrow(ci)){
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef[i], HR = hr[i], 
                                         CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
  }} else {
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef, HR = hr, 
                                         CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
  }

multi_os_cp <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1:5,11,6,7,10,8,9),]

res2$Indicator <- c(res2$Indicator[1:5],
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3", 
                    "N Stage\n(N1 Reference)", 
                    "N2", 
                    "N3")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.multivariate_cox_os_CP.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

# -----------多因素画图---------------
pdf(file = "07.multi_cox_DFS_cp.pdf", family = "Times", height = 10, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           # is.summary = c(TRUE,rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,2),TRUE,rep(FALSE,4)),
           is.summary = c(rep(FALSE,7),TRUE,rep(FALSE,2),TRUE,rep(FALSE,2)),
           # is.summary = c(rep(FALSE,16),TRUE,rep(FALSE,2),TRUE,rep(FALSE,2)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           # xticks = c(0, 1, 2, 4, 6, 8), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.4,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1, fontfamily = "Times"),
                          ticks=gpar(cex=0.8, fontface = "bold", fontfamily = "Times"),
                          xlab=gpar(cex = 1, fontface = "bold", fontfamily = "Times"),
                          title=gpar(cex = 1.25, fontface = "bold", fontfamily = "Times")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate",
           clip = c(0,7)) # 垂直于x轴的网格线，对应每个刻度
dev.off()
# -----------
# ----------DFS----------------
# -----------clinical----------------
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(clinical, collapse = "+")))
cox_more <- coxph(cox_data, data = data)

cox_zph <- cox.zph(cox_more)
names(cox_more$coefficients) <- clinical

mul_cox_result <- summary(cox_more)$coefficients

summary_cox <- summary(cox_more)

coef <- summary_cox$coefficients[, "coef"]
hr <- summary_cox$coefficients[, "exp(coef)"]
ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
p_value <- summary_cox$coefficients[, "Pr(>|z|)"]

results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

if (is.matrix(ci)){
  for (i in 1:nrow(ci)){
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef[i], HR = hr[i], 
                                         CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
  }} else {
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef, HR = hr, 
                                         CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
  }

multi_DFS_clinical <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res2 <- rbind(res) %>% as.data.frame()

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.multivariate_cox_dfs_clinical.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

# ------------------pathological-------------------
pathological <- c("pathological_grade","LVI","PNI","PFI","T_stage","N_stage")
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(pathological, collapse = "+")))
cox_more <- coxph(cox_data, data = data)

cox_zph <- cox.zph(cox_more)
names(cox_more$coefficients) <- c(pathological[1:4],
                                  "T2", 
                                  "T3",
                                  "N2", 
                                  "N3")

mul_cox_result <- summary(cox_more)$coefficients

summary_cox <- summary(cox_more)

coef <- summary_cox$coefficients[, "coef"]
hr <- summary_cox$coefficients[, "exp(coef)"]
ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
p_value <- summary_cox$coefficients[, "Pr(>|z|)"]

results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

if (is.matrix(ci)){
  for (i in 1:nrow(ci)){
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef[i], HR = hr[i], 
                                         CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
  }} else {
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef, HR = hr, 
                                         CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
  }

multi_DFS_pathological <- results$Variable[which(results$P_value < 0.05)]

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1:4,10,5,6,9,7,8),]

res2$Indicator <- c(res2$Indicator[1:4],
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3", 
                    "N Stage\n(N1 Reference)", 
                    "N2", 
                    "N3")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.multivariate_cox_dfs_pathological.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

# --------------cp-----------------
CP <- c("CA199","tumor_size","pathological_grade","LVI","PNI","PFI","T_stage","N_stage")
cox_data <- as.formula(paste0('Surv(DFS, recurrence)~', paste(CP, collapse = "+")))
cox_more <- coxph(cox_data, data = data)

cox_zph <- cox.zph(cox_more)
names(cox_more$coefficients) <- c(CP[1:6],
                                  "T2", 
                                  "T3",
                                  "N2", 
                                  "N3")

# cox_formula <- as.formula(paste("Surv(OS, censor)~",
#                                 paste(rownames(cox_table)[cox_table[,3]>0.05],
#                                       collapse = "+")))
# cox_more_2 <- coxph(cox_formula, data = data)
mul_cox_result <- summary(cox_more)$coefficients

summary_cox <- summary(cox_more)

coef <- summary_cox$coefficients[, "coef"]
hr <- summary_cox$coefficients[, "exp(coef)"]
ci <- summary_cox$conf.int[, c("lower .95", "upper .95")]
p_value <- summary_cox$coefficients[, "Pr(>|z|)"]

results <- data.frame(Variable = character(), Coefficient = numeric(), 
                      HR = numeric(), CI_lower = numeric(), CI_upper = numeric(), P_value = numeric(), 
                      stringsAsFactors = FALSE)

if (is.matrix(ci)){
  for (i in 1:nrow(ci)){
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef[i], HR = hr[i], 
                                         CI_lower = ci[,1][i], CI_upper = ci[,2][i], P_value = p_value[i]))
  }} else {
    results <- rbind(results, data.frame(Variable = rownames(summary_cox$coefficients)[i], Coefficient = coef, HR = hr, 
                                         CI_lower = ci[1], CI_upper = ci[2], P_value = p_value))
  }

res <- data.frame(p.value=results$P_value,
                  HR=results$HR,
                  HR.95L=results$CI_lower,
                  HR.95H=results$CI_upper,
                  Indicator=results$Variable)

res.ref <- c(1, 1, 1, 1, 1)

res2 <- rbind(res,res.ref,res.ref) %>% as.data.frame()
res2 <- res2[c(1:6,11,7,8,12,9,10),]

res2$Indicator <- c(res2$Indicator[1:6],
                    "T Stage\n(T1 Reference)", 
                    "T2", 
                    "T3", 
                    "N Stage\n(N1 Reference)", 
                    "N2", 
                    "N3")

rownames(res2) <- res2$Indicator
write.table(res2, file = "01.multivariate_cox_dfs_CP.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))

hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
# hz[c(4,8,11,14,18)] <- ""
tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.001,
                                      "< 0.001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

