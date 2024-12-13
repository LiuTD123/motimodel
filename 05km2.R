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

ROC_all$RS_group <- ifelse(ROC_all[,"RS"] > summary(cut_RS)[1,1],'High','Low' )
ROC_all$DL_group <- ifelse(ROC_all[,"DL"] > summary(cut_DL)[1,1],'High','Low' )
ROC_all$cp_group <- ifelse(ROC_all[,"cp"] > summary(cut_cp)[1,1],'High','Low' )
ROC_all$ajcc_group <- ifelse(ROC_all[,"ajcc"] > summary(cut_ajcc)[1,1],'High','Low' )

save(ROC_all,file = paste0("ROC_all_",datatype,"_riskgroup.RData"))

kmfit_rs <- survfit(Surv(futime, fustat) ~ RS_group, data = ROC_all)
kmfit_dl <- survfit(Surv(futime, fustat) ~ DL_group, data = ROC_all)
kmfit_cp <- survfit(Surv(futime, fustat) ~ cp_group, data = ROC_all)
kmfit_ajcc <- survfit(Surv(futime, fustat) ~ ajcc_group, data = ROC_all)

customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

# 定义参数列表
fit_list <- list(kmfit_rs, kmfit_cp, kmfit_dl, kmfit_ajcc)

palette_colors <- c("#f16c23","#2b6a99")
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
  train_km <- ggsurvplot(fit,
                         pval = TRUE, 
                         pval.method = T,
                         conf.int = F,
                         legend.labs=c("High risk","Low risk" ),
                         legend.title="Risk Score",
                         xlab = datatype,
                         title=paste0(str_to_title(dataset)," KM"),
                         font.main = c(15,"bold"),
                         risk.table = TRUE, 
                         risk.table.col = "strata", 
                         linetype = "strata", 
                         surv.median.line = "hv",
                         font.family = "Times",
                         risk.table.y.text.col = T,
                         risk.table.y.text = T,
                         risk.table.height = 0.35,
                         ggtheme = theme_bw(), 
                         palette = palette_colors)
  train_km$plot <- train_km$plot + labs(
    title    = paste0("Survival Curves in ",title),
    subtitle = "Based on Kaplan-Meier Estimates"
  )
  train_km$table <- train_km$table + labs(
    caption  = paste0("Created with ",str_to_title(dataset)," Data")
  )
  train_km <- customize_labels(
    train_km,
    font.title    = c(16, "bold"),
    font.subtitle = c(15, "bold.italic"),
    font.caption  = c(14, "plain", "orange"),
    font.x        = c(14, "bold.italic"),
    font.y        = c(14, "bold.italic"),
    font.xtickslab = c(12, "plain")
  )
  pvalue <- stringr::str_extract(train_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]],
                                 "\\d.*")
  train_km[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(pvalue)))
  
  pdf_file_name <- paste0(title, ".pdf")
  pdf(pdf_file_name, family = "Times", height = 6, width = 6, onefile = F)
  print(train_km)
  dev.off()
}

# -------------------dis--------------
varlist <-c("RS","DL","cp","ajcc")

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