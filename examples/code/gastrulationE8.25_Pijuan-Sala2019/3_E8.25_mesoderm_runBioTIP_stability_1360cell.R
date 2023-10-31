# This code doing:
# extract the QUANTC-assigned 4 clusters for 1362 cells
# run BioTIP on these 1362 cells
# plot resutls
# last update 1/14/2022
# by Holly yang

setwd('F:/projects/scRNA/results/GSE87038_gastrulation/uncorrected/E8.25_mesoderm_robustness')
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/stability_score.R') 

# the following funciton has been updated on github, but before updating the package to Bioconductor, you need to call it locally to test
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/ForJennifer/optimize.sd_selection.R') 
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/BioTIP.wrap.R')

# normalized data were downloaded and reanalyzed
# refer to GSE87038_E8.25_normalize.R
library(BioTIP)
library(dplyr)
library(SingleCellExperiment)

library(ggplot2)
library(gridExtra)  # required by grid.arrange()
library(ggpubr)  # required by p-values on bar



###################################################
## 1) load the R object ##
## refer to 2_QuanTC_simulation_cluster_notused.R
## focusing on 1.3k cells analyzed by QuanTC
## refer to BioTIP_E8.25_mesoderm_add_QuanTC.clusters.R
###################################################

load('sce_E8.25_HEP.RData')
sce
# sceclass: SingleCellExperiment 
# dim: 3073 1362 
# metadata(0):
#   assays(2): counts logcounts
# rownames(3073): Phlda2 Myl7 ... Hmcn1 Tfdp2
# rowData names(2): ENSEMBL SYMBOL
# colnames(1362): cell_63240 cell_63242 ... cell_95715 cell_95717
# colData names(21): cell barcode ... C_Leiden_0.8 C_Leiden_1.2
# reducedDimNames(5): pca.corrected umap TSNE UMAP force
# altExpNames(0):


load(file='CTS_ShrinkM_E8.25_HEP.RData')
dim(M) #3073 3073
all(rownames(sce) == rownames(M)) #TRUE
  
table(colData(sce)$label)
# 3  13  10   6  15   7 
# 381 223 278 283  60 137 
#colData(sce)$label <- factor(sce$label, levels=c('3', '13', '10','6','15', '7'))
class(colData(sce)$label)
#[1] "factor"



######### 1) setting parameters,   ######################
localHVG.preselect.cut = 0.1 # A positive numeric value setting how many propotion of global HVG to be selected per cell cluster
getNetwork.cut.fdr = 0.2  # to construct RW network and extract co-expressed gene moduels
getTopMCI.gene.minsize = c(10,20,30,40)  # min number of genes in an identified CTS
getTopMCI.n.states = 5  # A number setting the number of states to check the significance of DNB score (i.e. MCI) 
                        # This parameter can be enlarge when inputting more clusters
MCIbottom = 2:4  # A number setting the threshold to prioritize the initially selected CTS candidates.
# In our experiment, a number between 2 and 4 generated expected resutls.

############### run BioTIP  ##################
n.scability <- 20
downsize <- 0.95  
n1 <- nrow(sce)
n2 <- ncol(sce)

for(m in 1:length(getTopMCI.gene.minsize)){
  for(i in 1:length(MCIbottom)){
    set.seed(102020)
    ## get n.scability sets of predictions
    res <- list()  
    flag = counts = 0
    repeat{ 
      cat(counts, '\t')
      counts = counts+1 # control the maxmum trials for searching two sets of signficant CTSs, given the current parameter settings
      sce.dnsize <- sce[sample(n1, floor(n1*downsize)), sample(n2, floor(n2*downsize))]
      samplesL <- split(rownames(colData(sce.dnsize)), f = sce.dnsize$label)
      this.run <- BioTIP.wrap(sce.dnsize, samplesL, subDir='stability', 
                              getTopMCI.n.states=getTopMCI.n.states, 
                              getTopMCI.gene.minsize= floor(getTopMCI.gene.minsize[m] * downsize), 
                              MCIbottom=MCIbottom[i],
                              localHVG.preselect.cut=localHVG.preselect.cut, 
                              getNetwork.cut.fdr=getNetwork.cut.fdr, 
                              M=M,
                              verbose=FALSE, plot=FALSE)   
      if(!is.null(this.run) ) { 
        if(length(which(this.run$significant))>0)  { # control at least 1 significant CT to be identified !!!
          flag = flag +1
          res[[flag]] <- this.run
          #rm(this.run)
        }
      }
      if(flag ==n.scability | counts == 100) {
        break
      } 
    } # end repeat
    
    counts-flag  # 1, the number of runs without significant  
    save(res, flag, counts, file=paste0('stability_1360cell/BioTIP_res_20runs_MCIbottom',MCIbottom[i],
                                        '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    
  }
}

###########################################################
## evaluate the stability of CT detection
## using the robustly detected C13, C15, C6, and C7 as reference
###########################################################
df.CT =NULL
for(m in 1:length(getTopMCI.gene.minsize)){
  CT.list <- list()
  Freq <- matrix(0, nrow=nlevels(sce$label), ncol=length(MCIbottom))#!!!
  rownames(Freq) <- levels(sce$label) #!!!
  colnames(Freq) <- MCIbottom 
  for(i in 1:length(MCIbottom)){
    load(file=paste0('stability_1360cell/BioTIP_res_20runs_MCIbottom',MCIbottom[i],
                     '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    ## count which cluster has been identified repeatily ###########
    CT <- list()
   # for(j in 1:n.stability)  
    for(j in 1:length(res)){ # in case n.stability>length(res), i.e. some run has no isgnificant prediction
      CT[[j]] <- res[[j]]$significant[which(res[[j]]$significant)] %>% names() %>% unique()
    }
    tmp <- table(unlist(CT))
    Freq[names(tmp),i] <- tmp
    #Freq <- cbind(Freq, tmp)
    CT.list[[i]] <- CT
    } 

 # Freq <- Freq/length(res)
  names(CT.list) <- MCIbottom 
  
  df <- data.frame(Frequency = as.numeric(Freq),
                   clusterID = rep(rownames(Freq), ncol(Freq)),
                   MCIbottom = lapply(as.character(MCIbottom ), function(x) rep(x,nrow(Freq))) %>% unlist()
  )
  head(df)
  #   Frequency clusterID MCIbottom
  # 1    10       C13         2
  # 2     1        C3         2
  # 3     5        C5         2
  # 4     2        C7         2
  # 5     3        C8         2
  # 6    12        C6         2
  df$clusterID <- factor(df$clusterID, levels =levels(sce$label)) #!!!
  df$moduleSize <- rep(getTopMCI.gene.minsize[m], nrow(df))
  
  df.CT <- rbind(df.CT, df)  
 # save(CT.list, file= paste0('stability_1360cell/CT.detection_Modulesize',getTopMCI.gene.minsize[m],'.Rdata'))
}
df.CT$Frequency = df.CT$Frequency/n.scability
save(df.CT, file= 'stability_1360cell/CT.detection_variable.MCIbottom_Modulesize.Rdata')


#############################################################
## plot for variable Modulesize, when MCIbottom==2         ##
#############################################################
load(file= 'stability_1360cell/CT.detection_variable.MCIbottom_Modulesize.Rdata')
df.CT$moduleSize <- as.character(df.CT$moduleSize )
 ggplot(data=subset(df.CT, MCIbottom==2),
        aes(x=clusterID, y=Frequency, group=moduleSize)) +
    geom_line(aes(color=moduleSize), size=1.5)+
    geom_point(aes(color=moduleSize)) 
dev.copy2pdf(file='stability_1360cell/CT.detection_variable.MCIbottom2_Modulesize.pdf')
 



###########################################################
## evaluate the stability of CTS identification
## for the  C13  when MCIbottom==2, we did two estimations:
## method 1) we generate a proxy 'golden standard' of CTS  -- used                                               
## which were 16 genes identifed by at least 3 out of eight identfications for early HEP                             
## method 2) we used the original CTS_C13 
## refer to 2.3_E8.25_mesoderm_quantify_robustness.R
###########################################################
MCIbottom
#[1] 2 3 4
i=1  # simpliy check when MCIbottom==2 

GS.CTS <- read.table('GS.CTS.eHEP_at.least.3identificaitons_fr.8predicitons.txt')
GS.CTS2 <- read.table('GS.CTS.lHEP_at.least.4identificaitons_fr.8predicitons.txt')
GS.CTS3 <- read.table('GS.CTS.EP_at.least.4identificaitons_fr.7predicitons.txt')
GS.CTS4 <- read.table('GS.CTS.HP_at.least.5identificaitons_fr.9predicitons.txt')
ref <- list('13'=GS.CTS[,1], '15'=GS.CTS2[,1], '6'=GS.CTS3[,1], '7'=GS.CTS4[,1] )
rm(GS.CTS, GS.CTS2, GS.CTS3,GS.CTS4)

getTopMCI.gene.minsize <- getTopMCI.gene.minsize[1:3] # simplify the dicussion for only 3 options

for(k in 1:length(ref)){
  CT.ref <- ref[k]

  F1.CT <- F1.CT.ctl <- N.module <- list()
  for(m in 1:length(getTopMCI.gene.minsize)){
    load(file=paste0('stability_1360cell/BioTIP_res_20runs_MCIbottom',MCIbottom[i],
                     '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    F1.CTS.CT = F1.bk.CT = n.module <- array()
    for(j in 1:length(res)){
      # evaluate the F1 score for two CTS identificitons
      set2 <- res[[j]]$CTS.candidate[which(res[[j]]$significant)] 
      F1.CTS.CT[j] <- F_score.CTS(CT.ref, set2, weight=TRUE)$F1
      # negative control
      n <- lengths(set2)
      random.set2 <- lapply(n, function(x) sample(rownames(sce), x))
      F1.bk.CT[j] <- F_score.CTS(CT.ref, random.set2, weight=TRUE)$F1
      n.module[j] <- n[names(CT.ref)]
    } 
    F1.CT[[m]] =  F1.CTS.CT
    F1.CT.ctl[[m]] =  F1.bk.CT
    N.module[[m]] = n.module
  }
  names(F1.CT) <- names(F1.CT.ctl) <- names(N.module) <- getTopMCI.gene.minsize
  
  save(F1.CT, F1.CT.ctl, N.module,
       file=paste0('stability_1360cell/F1.CTS_MCIbottom2_GS.C',names(ref)[k],'.RData'))
  }

# plot 
{
  load('stability_1360cell/F1.CTS_MCIbottom2_GS.C13.RData')
  F1.C13 <- F1.CT 
  F1.C13.ctl <- F1.CT.ctl
  lapply(N.module, function(x) mean(x, na.rm=T)) %>% unlist()
  # 10       20       30 
  # 28.70000 35.15789 36.00000 
  load('stability_1360cell/F1.CTS_MCIbottom2_GS.C15.RData')
  F1.C15 <- F1.CT 
  F1.C15.ctl <- F1.CT.ctl
  lapply(N.module, function(x) mean(x, na.rm=T)) %>% unlist()
  # 10       20       30 
  # 15.25000 25.63636 30.25000 
  load('stability_1360cell/F1.CTS_MCIbottom2_GS.C6.RData')
  F1.C6 <- F1.CT 
  F1.C6.ctl <- F1.CT.ctl
  lapply(N.module, function(x) mean(x, na.rm=T)) %>% unlist()
  # 10       20       30 
  # 59.83333 63.31579 72.55000  >20% of the 3k global HVG
  load('stability_1360cell/F1.CTS_MCIbottom2_GS.C7.RData')
  F1.C7 <- F1.CT 
  F1.C7.ctl <- F1.CT.ctl
  lapply(N.module, function(x) mean(x, na.rm=T)) %>% unlist()
  # 10       20       30 
  # 17.90000 27.60000 34.09091   
  
  lapply(F1.C13, mean) %>% unlist()
  #        10         20         30 
  # 0.06079849 0.06412630 0.04162247 

  df <- data.frame(F1 = c(unlist(F1.C13.ctl),unlist(F1.C13),
                          unlist(F1.C15.ctl),unlist(F1.C15),
                          unlist(F1.C6.ctl),unlist(F1.C6),
                          unlist(F1.C7.ctl),unlist(F1.C7) ),
                   minModuleSize = c(rep(getTopMCI.gene.minsize, lengths(F1.C13.ctl)),
                                     rep(getTopMCI.gene.minsize, lengths(F1.C13)),
                                     rep(getTopMCI.gene.minsize, lengths(F1.C15.ctl)),
                                     rep(getTopMCI.gene.minsize, lengths(F1.C15)),
                                     rep(getTopMCI.gene.minsize, lengths(F1.C6.ctl)),
                                     rep(getTopMCI.gene.minsize, lengths(F1.C6)),
                                     rep(getTopMCI.gene.minsize, lengths(F1.C7.ctl)),
                                     rep(getTopMCI.gene.minsize, lengths(F1.C7)) ),
                   CT = c(rep('C13', sum(lengths(F1.C13.ctl))+sum(lengths(F1.C13))),
                          rep('C15', sum(lengths(F1.C15.ctl))+sum(lengths(F1.C15))),
                          rep('C6', sum(lengths(F1.C6.ctl))+sum(lengths(F1.C6))),
                          rep('C7', sum(lengths(F1.C7.ctl))+sum(lengths(F1.C7))) ),
                   CTS = c(rep('ctl', sum(lengths(F1.C13.ctl))), rep('CTS', sum(lengths(F1.C13))),
                           rep('ctl', sum(lengths(F1.C15.ctl))), rep('CTS', sum(lengths(F1.C15))),
                           rep('ctl', sum(lengths(F1.C6.ctl))), rep('CTS', sum(lengths(F1.C6))),
                           rep('ctl', sum(lengths(F1.C7.ctl))), rep('CTS', sum(lengths(F1.C7))) ) 
  )
  df$minModuleSize <- as.character(df$minModuleSize) 
  df$CTS <- factor(df$CTS, levels = c('ctl','CTS')) 
  df$color.lab <- paste0(df$minModuleSize,df$CTS)
  df$color.lab <- factor(df$color.lab, levels=c('10ctl','10CTS','20ctl','20CTS','30ctl','30CTS')) #, '40ctl','40CTS'))
  df$CT <- factor(df$CT, levels =c('C13','C15','C6','C7'))
  
  df$Norm.F1 = Normalize.F1(df$F1)
  
  # To ensure that easy and difficult CT state have equal influence on the final score, 
  # we normalized the scores at each CT state across the different MinModuleSize setting.
  
  library(ggbeeswarm)
  library(grid)
  library(gridExtra)
 
  # plot w C6, used 
  {
  p1 <- ggplot(data=df, aes(x=CT, y=F1, color=color.lab)) +
    geom_violin(position = position_dodge(width = 0.5)) +
    geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
    xlab('predicted for cluster') + ylab('F1 score for the CTS') +
    scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue", "grey","red")) 
  
  p2 <- ggplot(data=df, aes(x=CT, y=Norm.F1, color=color.lab)) +
    geom_violin(position = position_dodge(width = 0.5)) +
    geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
    xlab('predicted for cluster') + ylab('Normalized F1 score for the CTS') +
    scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue", "grey","red")) 
  
  pdf(file=paste0('stability_1360cell/CTS.F1_variable.ModuleSize_GS_4CTs.pdf'))
  gridExtra::grid.arrange(
    p1, p2, top = "CTS stability", bottom = "different ModuleSize settings"
    ,ncol=1)
  dev.off()
}
} 
dim(sce)  #3073
## prsent the statistics using proxy gold standard here an in the figure
# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

tmp <- compare_means(Norm.F1 ~ color.lab, data=df, group='CT', method='t.test')
subset(tmp, CT=='C13' & group1=='10ctl' & group2=='10CTS')
#   1  C13   Norm.F1 10ctl  10CTS  5.66e-8 2.3e-6 5.7e-08  **** 
subset(tmp, CT=='C13' & group1=='20ctl' & group2=='20CTS')
#   1 C13   Norm.F1 20ctl  20CTS  9.63e-6 3.8e-4 9.6e-06  ****  
subset(tmp, CT=='C13' & group1=='30ctl' & group2=='30CTS')
# 1 C13   Norm.F1 30ctl  30CTS  0.00614  0.16 0.00614  **  

subset(tmp, CT=='C15' & group1=='10ctl' & group2=='10CTS')
# <chr> <chr>   <chr>  <chr>     <dbl>   <dbl> <chr>    <chr>    <chr> 
#   1  C15   Norm.F1 10ctl  10CTS  3.26e-8 1.4e-6 3.3e-08  **** 
subset(tmp, CT=='C15' & group1=='20ctl' & group2=='20CTS')
#   1 C15   Norm.F1 20ctl  20CTS  0.00134 0.041 0.00134  **  
subset(tmp, CT=='C15' & group1=='30ctl' & group2=='30CTS')
# 1C15   Norm.F1 30ctl  30CTS  0.0930     1 0.09296  ns    


subset(tmp, CT=='C6' & group1=='10ctl' & group2=='10CTS')
#   1 C6    Norm.F1 10ctl  10CTS  1.08e-10 0.0000000051 1.1e-10  ****     T-test
subset(tmp, CT=='C6' & group1=='20ctl' & group2=='20CTS')
#   1 C6    Norm.F1 20ctl  20CTS  8.93e-14 4.8e-12 8.9e-14  ****     T-test
subset(tmp, CT=='C6' & group1=='30ctl' & group2=='30CTS')
#   1 C6    Norm.F1 30ctl  30CTS  5.53e-37 3.2e-35 < 2e-16  ****     T-test

subset(tmp, CT=='C7' & group1=='10ctl' & group2=='10CTS')
#   1 C7    Norm.F1 10ctl  10CTS  1.67e-12 8.3e-11 1.7e-12  ****     T-test
subset(tmp, CT=='C7' & group1=='20ctl' & group2=='20CTS')
#   1 C7    Norm.F1 20ctl  20CTS  3.37e-17 1.9e-15 < 2e-16  ****     T-test
subset(tmp, CT=='C7' & group1=='30ctl' & group2=='30CTS')
#   1 C7    Norm.F1 30ctl  30CTS  0.000241 0.0084 0.00024  ***      T-test

