rm(list = ls())
setwd("E:/workdir/08multimodel/results/DFS")
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

if (!requireNamespace("fastAdaboost", quietly = TRUE))
  devtools::install_github("souravc83/fastAdaboost")

# if (!requireNamespace("Mime", quietly = TRUE))
#   devtools::install_github("l-magnificence/Mime")

library(Mime1)
library(readxl)

# -------OS----------
# ----------------训练集--------------------
data_train <- read_excel("../../data_train2.xlsx")
# features <- read.csv("../../03cox/cox_sigfeature.csv",header = F)
features <- c("CA199","tumorsize","pathologicalgrade","PNI","Tstage","Nstage")

# data_train <- data_train[,c("id","OS","censor","OS_risk",features)]
data_train <- data_train[,c("id","DFS","recurrence","DFS_risk",features)]
colnames(data_train)[1:4] <- c("ID","OS.time","OS","DL")

# fatorfeatures <- c("CA199","pathologicalgrade","PNI","Tstage","Nstage")
# for (var in fatorfeatures){
#   data_train[[var]] <- factor(data_train[[var]])
# }

# ----------------验证集--------------------
data_test <- read_excel("E:/workdir/08multimodel/data_test2.xlsx")
# colnames(data_test)[c(9,8,10,6)] <- c("OS","OS.time","ID","DL")
colnames(data_test)[c(5,4,10,2)] <- c("OS","OS.time","ID","DL")

data_test <- data_test[,colnames(data_train)]

# for (var in fatorfeatures){
#   data_test[[var]] <- factor(data_test[[var]])
# }
# ---------外部--------------
data_test_e <- read_excel("E:/workdir/08multimodel/data_all_external2.xlsx")
# colnames(data_test_e)[c(9,8,10,6)] <- c("OS","OS.time","ID","DL")
colnames(data_test_e)[c(5,4,10,2)] <- c("OS","OS.time","ID","DL")

data_test_e <- data_test_e[,colnames(data_train)]
# for (var in fatorfeatures){
#   data_test_e[[var]] <- factor(data_test_e[[var]])
# }
# -------group# ------------------------------------------
# load("./exampledata/Example.cohort.Rdata") # 生存数据与基因表达信息
# load("./exampledata/genelist.Rdata")
# list_train_vali_Data[[5]][1:5,1:5]$OS
#                 ID    OS.time OS   MT-CO1   MT-CO3
#60  TCGA.DH.A66B.01 1281.65322  0 13.77340 13.67931
#234 TCGA.HT.7607.01   96.19915  1 14.96535 14.31857
#42  TCGA.DB.A64Q.01  182.37755  0 13.90659 13.65321
#126 TCGA.DU.8167.01  471.97707  0 14.90695 14.59776
#237 TCGA.HT.7610.01 1709.53901  0 15.22784 14.62756
# 其中list_train_vali_Data是含有2个数据集的列表，每个数据集的第一列为ID ，2-3列为生存信息（OS.time ，OS） ，后面为基因表达量。

list_train_vali_Data <- list("Train" = data_train,
                             "Test" = data_test,
                             "External" = data_test_e
                             )

save(list_train_vali_Data, file = "train_vali_data.RData")

load("train_vali_data.RData")
# 二 构建预后模型
# 1. 构建101机器学习模型组合
genelist <- colnames(data_train)[4:10]
colnames(data_train)

res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Train,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.5,
                       candidate_genes = genelist,
                       mode = 'all',
                       nodesize =5,
                       seed = 599999999 )

save(res, file = "res_all_ML.Dev.DFS.RData")

load("res_all_ML.Dev.DFS.RData")

# 提取系数
ml <- res$ml.res
ml <- ml$`StepCox[forward] + RSF`
modelgenes <- ml$`StepCox[both] + GBM`$fit$var.names

save(ml,file = "ml_details.RData")
write.csv(modelgenes,file = "modelgenes.csv")

# ML.Dev.Prog.Sig() 可选 all, single 和 double三种模式. all 为所有10种算法 以及 组合 . 
# single 为用10种算法中的一种.
# double 为两种算法的组合，一般情况下使用 all 模式.
# 默认情况下 unicox.filter.for.candi 为 T , 会先对训练集进行单因素cox分析，unicox_p_cutoff 显著的基因会用于构建预后模型.

# 如果使用自己数据的时候，需要注意：
# （1）替换自己数据注意前三列的要求，且将多个数据集以列表形式存储
# （2）分析之前最好先确认 所有数据集中是否 有基因集列表中的所有基因 ，减少报错。
# （3）种子数确定好，会有一些小的影响 。

data2 <- data2 %>% 
  dplyr::select(ID , OS.time , OS, genelist)
# 通过View(ML.Dev.Prog.Sig) 检查函数，设置unicox.filter.for.candi = T 后会先做单因素cox分析

# 2. C-index 展示
# 示例数据list_train_vali_Data 为2个数据集的list，结果图中队列为2个，最后两列为Cindex的均值，这也就是机器学习模型组合文献中的主图。
pdf("cindex_dis.pdf",w=9,h=18)
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[1],
               order = names(list_train_vali_Data),
               width = 0.35
)
dev.off()

cindex <- res$Cindex.res
write.csv(cindex,file = "cindex.csv",quote = F)

# 3. 查看指定模型的结果
# 假设我们选择第一个模型（StepCox[forward] + plsRcox） ，可以单独查看该模型下各个数据集的cindex表现

method_select = "SuperPC"

pdf("cindex_dis_select.pdf",w=5,h=5)
cindex_dis_select(res,
                  model= method_select,
                  order= names(list_train_vali_Data))
dev.off()

# 也可以查看该模型下各个数据集的KM曲线情况
survplot <- vector("list",3) 
pdf("km_.pdf",w=8,h=7)
for (i in c(1:3)) {
  survplot[[i]]<-rs_sur(res, model_name = method_select,
                        dataset = names(list_train_vali_Data)[i],
                        # color=c("blue","green"),
                        median.line = "hv",
                        # cutoff = 0.5,
                        conf.int = T,
                        xlab="Day",
                        pval.coord=c(1000,0.9)
                              )
  
}
dev.off()

# rs_sur(res, model_name = method_select,
#        dataset = names(list_train_vali_Data)[6],
#        # color=c("blue","green"),
#        median.line = "hv",
#        cutoff = 0.5,
#        conf.int = T,
#        xlab="Day",
#        pval.coord=c(1000,0.9)
# )

pdf("km_combine.pdf",w=16,h=5)
aplot::plot_list(gglist=survplot,ncol=3)
dev.off()

# 提取模型RS结果
# 这里有个很重要的点是要提取指定模型下的RS结果，然后就可以根据自己的需求重新绘制KM 以及 独立预后分析，森林图，列线图等其他分析了。
# 结果都在res中，根据str(res)知道对应的信息，提取即可

head(res$riskscore$`StepCox[forward] + RSF`[[1]])
head(res$riskscore$`StepCox[forward] + RSF`[[2]])

trainriskscore <- res$riskscore[[method_select]][[1]]
testriskscore <- res$riskscore[[method_select]][[2]]
externalriskscore <- res$riskscore[[method_select]][[3]]

save(trainriskscore,testriskscore,externalriskscore,file = "riskscore_superpc.RData")
# 提取特定模型的基因
ml <- res$ml.res
modelgenes <- ml$`StepCox[both] + GBM`$terms
genes <- colnames(modelgenes)

save(modelgenes,genes,file = "modelgenes.RData")

# 4. AUC结果
# 计算每个模型的1年，3年，5年 的 auc值 ，并可视化所有模型的1年auc结果
cutoff_1 <- 1
cutoff_2 <- 2
cutoff_3 <- 3
cutoff_4 <- 4
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Train"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = "all",
                             AUC_time = cutoff_1,
                             auc_cal_method="KM")
all.auc.2y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Train"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = "all",
                             AUC_time = cutoff_2,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Train"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = "all",
                             AUC_time = cutoff_3,
                             auc_cal_method="KM")
all.auc.4y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Train"]],
                             inputmatrix.list = list_train_vali_Data,
                             mode = "all",
                             AUC_time = cutoff_4,
                             auc_cal_method="KM")
save.image("auc_over.RData")

pdf(file = paste0("roc-",cutoff_1,"year.pdf"),w=9,h=18)
auc_dis_all(all.auc.1y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off()

pdf(paste0("roc-",cutoff_2,"year.pdf"),w=9,h=18)
auc_dis_all(all.auc.2y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off()

pdf(paste0("roc-",cutoff_3,"year.pdf"),w=9,h=18)
auc_dis_all(all.auc.3y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off()

pdf(paste0("roc-",cutoff_4,"year.pdf"),w=9,h=18)
auc_dis_all(all.auc.4y,
            dataset = names(list_train_vali_Data),
            validate_set=names(list_train_vali_Data)[-1],
            order= names(list_train_vali_Data),
            width = 0.35,
            year=1)
dev.off()

# 绘制选定模型下的auc曲线
roc_vis(all.auc.5y,
        model_name = method_select,
        dataset = names(list_train_vali_Data),
        order= names(list_train_vali_Data),
        anno_position=c(0.65,0.55),
        year=1)

pdf("auc_histon.pdf",w=7,h=5)
auc_dis_select(list(all.auc.1y,all.auc.3y,all.auc.5y),
               model_name= method_select,
               dataset = names(list_train_vali_Data),
               order= names(list_train_vali_Data),
               year=c(cutoff_1,cutoff_2,cutoff_3))
dev.off()

# 5. 模型比较
# 该包还提供了和之前文献报道的预后模型比较的函数，当然只提供了胶质瘤的。
# 那如果你做的是其他癌种呢？可以通过查看函数了解是怎样的输入形式，然后就做对应的替换后就可以分析

cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = F,
                                             type.sig = c('Glioma','LGG','GBM'),
                                             list_input_data = list_train_vali_Data)

pdf("model_compare_dataset5_rsf.pdf",w=9,h=16)
cindex_comp(cc.glioma.lgg.gbm,
            res,
            model_name="StepCox[forward] + RSF",
            dataset=names(list_train_vali_Data))
dev.off()
