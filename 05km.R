rm(list = ls())

dataset = "external"
datatype = "OS"
foldpath <- paste0("E:/workdir/08multimodel/05km/",datatype,"/",dataset)

if (!dir.exists(foldpath)){
  dir.create(foldpath, recursive = TRUE)
}
setwd(foldpath)
library(survival)
library(survminer)
load(paste0("../../../04auccompare/COX/",datatype,"/",dataset,"ROC_all.RData"))

# colnames(ROC_all) <- c("ID","RS","DL","fustat","futime","cp","ajcc")

cut_RS <- survminer::surv_cutpoint(ROC_all, #数据集
                                   minprop = 0.25,
                                   time = "futime", #生存时间
                                   event = "fustat", #生存状态
                                   variables = "RS"  #需要计算的数据列名
)
cut_RS

cut_DL <- survminer::surv_cutpoint(ROC_all, #数据集
                                   minprop = 0.25,
                                   time = "futime", #生存时间
                                   event = "fustat", #生存状态
                                   variables = "DL"  #需要计算的数据列名
)
cut_DL

cut_cp <- survminer::surv_cutpoint(ROC_all, #数据集
                                   minprop = 0.25,
                                   time = "futime", #生存时间
                                   event = "fustat", #生存状态
                                   variables = "cp"  #需要计算的数据列名
)
cut_cp

cut_ajcc <- survminer::surv_cutpoint(ROC_all, #数据集
                                   minprop = 0.25,
                                   time = "futime", #生存时间
                                   event = "fustat", #生存状态
                                   variables = "ajcc"  #需要计算的数据列名
)
cut_ajcc

ROC_all$RS_group <- as.data.frame(ifelse(ROC_all[,"RS"] > summary(cut_RS)[1,1],'High','Low' ))[,1]
# ROC_all$RS_group <- factor(ROC_all$RS_group, 
#                            levels = c("High", "Low"), 
#                            labels = c("DHR", "DLR"))
ROC_all$DL_group <- as.data.frame(ifelse(ROC_all[,"DL"] > summary(cut_RS)[1,1],'High','Low' ))[,1]
# ROC_all$DL_group <- factor(ROC_all$DL_group, 
#                            levels = c("High", "Low"), 
#                            labels = c("DHR", "DLR"))
ROC_all$cp_group <- as.data.frame(ifelse(ROC_all[,"cp"] > summary(cut_RS)[1,1],'High','Low' ))[,1]
# ROC_all$cp_group <- factor(ROC_all$cp_group, 
#                            levels = c("High", "Low"), 
#                            labels = c("DHR", "DLR"))
ROC_all$ajcc_group <- as.data.frame(ifelse(ROC_all[,"ajcc"] > summary(cut_RS)[1,1],'High','Low' ))[,1]
# ROC_all$ajcc_group <- factor(ROC_all$ajcc_group, 
#                            levels = c("High", "Low"), 
#                            labels = c("DHR", "DLR"))

save(ROC_all,file = paste0("ROC_all_",datatype,"_riskgroup.RData"))

palette_colors <- c("#f16c23","#2b6a99")

kmfit_rs <- survfit(Surv(futime, fustat) ~ RS_group, data = ROC_all)
kmfit_dl <- survfit(Surv(futime, fustat) ~ DL_group, data = ROC_all)
kmfit_cp <- survfit(Surv(futime, fustat) ~ cp_group, data = ROC_all)
kmfit_ajcc <- survfit(Surv(futime, fustat) ~ ajcc_group, data = ROC_all)

custom_theme <- function() {
  theme_survminer(
    font.risk.table = c(10, "bold"),
    risk.table.y.text.col = TRUE, # 分组名称颜色匹配
    font.tickslab = c(10, "bold"), # 轴文字加粗
    font.title = c(12, "bold"),
    font.subtitle = c(10, "italic"),
    font.legend = c(10, "bold")
  )
}

# 定义参数列表
fit_list <- list(kmfit_rs, kmfit_cp, kmfit_dl, kmfit_ajcc)
i = 1
for (i in 1:length(fit_list)) {
  fit <- fit_list[[i]]
  # 根据拟合对象的名称设置标题
  if (i == 1) {
    title <- "Merge KM"
  } else if (i == 2) {
    title <- "CP KM"
  } else if (i == 3) {
    title <- "DL KM"
  } else {
    title = "AJCC KM"
  }
  train_km <- ggsurvplot(
    fit, 
    # data = ROC_all,
    pval = TRUE, 
    risk.table = TRUE, 
    surv.median.line = "hv",
    palette = palette_colors, 
    legend.title = "Risk Group",
    title = "Derivation", 
    ylab = paste0("Cumulative ",datatype), 
    xlab = "Time (Days)", 
    break.x.by = 365,
    risk.table.height = 0.3,
    risk.table.y.text = TRUE,  
    risk.table.y.text.col = TRUE, # 启用分组名称颜色匹配
    theme = custom_theme(),
    risk.table.col = "strata"
  )

  pdf_file_name <- paste0(title, ".pdf")
  pdf(pdf_file_name, family = "Times", height = 6, width = 6, onefile = F)
  print(train_km)
  dev.off()
}

# -------------------dis--------------
varlist <-c("RS","DL","cp","ajcc")
var <- "cp"
for (var in varlist){
  vargroup <- paste0(var,"_group")
  title <- var
  res <- survminer::surv_cutpoint(ROC_all, #数据集
                                     minprop = 0.25,
                                     time = "futime", #生存时间
                                     event = "fustat", #生存状态
                                     variables = var  #需要计算的数据列名
  )
  risk_dis <- ggplot(ROC_all, aes(x=reorder(ROC_all$id, ROC_all[[var]]), y=ROC_all[[var]], color = ROC_all[[vargroup]])) +
    geom_point() +
    scale_color_manual(values = palette_colors) + 
    scale_x_discrete(breaks = rownames(ROC_all$id)[order(ROC_all[[var]])][c(1,50,100,150,200,250,300,350,400,450)],
                     labels = c(1,50,100,150,200,250,300,350,400,450),
                     expand = c(0.02,0)) +
    geom_vline(xintercept = nrow(ROC_all[which(ROC_all[[vargroup]]=="Low"),]) + 0.5, lty = 2) +
    # geom_hline(yintercept = summary(res)[1,1], lty =2) +
    labs(x = "Patients(increasing risk score)",
         y = "Risk Score",
         title = "Risk Score Distribution") + 
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1),
          legend.position = c(0,1),
          # legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(color = "black", size = .3),
          panel.grid = element_blank(),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 15),
          text = element_text(family = "Times", face = "bold"))
  ggsave(filename = paste0(var,"_riskScore_dis.pdf"), height = 3, width = 5, risk_dis)
}

for (var in varlist){
  vargroup <- paste0(var,"_group")
  title <- var
  res <- survminer::surv_cutpoint(ROC_all, #数据集
                                  minprop = 0.25,
                                  time = "futime", #生存时间
                                  event = "fustat", #生存状态
                                  variables = var  #需要计算的数据列名
  )
  if (datatype == "OS"){
    lable <- c("Dead","Alive")
  } else {
    lable <- c("Reccurence","DiseaseFree")
  }
  surv_stat <- ggplot(ROC_all, aes(x=reorder(ROC_all$id, ROC_all[[var]]),
                                              y=futime,
                                              color = factor(fustat,
                                                             levels = c(0,1),
                                                             labels = lable))) +
    geom_point() +
    scale_color_manual(values = palette_colors) +
    scale_x_discrete(breaks = ROC_all$id[order(ROC_all[[var]])][c(1,50,100,150,200,250,300,350,400,450)],
                     labels = c(1,50,100,150,200,250,300,350,400,450),
                     expand = c(0.02,0)) +
    geom_vline(xintercept = nrow(ROC_all[which(ROC_all[[vargroup]]=="Low"),]) + 0.5, lty = 2) +
    labs(x = "Patients(increasing risk score)",
         y = paste0(datatype," (days)"),
         title = "Distribution") + 
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.key = element_blank(),
          legend.justification = c(0,1),
          legend.position = c(0,1),
          # legend.margin = margin(c(-5,4,4,3)),
          legend.background = element_rect(color = "black", size = .3),
          panel.grid = element_blank(),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 15),
          text = element_text(family = "Times", face = "bold"))
  ggsave(filename = paste0(var,"_dis.pdf"), height = 3, width = 5, surv_stat)
}
