# This code doing:
# 1) Run BioTIP iterativly on 95% downsized cells and 95% dowmsized genes, 10 times
#
# 2) FOr each each pair of subsequent outputs, calcualte two matrics:
## Frequency of recapturing the original CT detections; normalized F1 scores for CTS identifications
#
#
# last update 1/28/2022
# by Holly yang

setwd('F:/projects/BioTIP/result/Bargaje2017_EOMES/hESC_Bargaje2017_robustness/')
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
## load the R object ##
## focusing on 1.3k cells of the HEP processes
## refer to BioTIP_E8.25_mesoderm_add.clusters.R
###################################################

load('sce_hESC.RData')
sce
# class: SingleCellExperiment 
# dim: 96 929
 
samplesL <- split(rownames(colData(sce)), f = colData(sce)$Consensus.Cluster)
(tmp = lengths(samplesL))
tmp <- tmp[which(tmp > 20)]

sce$Consensus.Cluster <- factor(as.vector(sce$Consensus.Cluster), 
                                        levels=c( "1" ,"10" ,"2",  "3",  "4",  "5",  "7" , "8",  "9"))

# logmat <- as.matrix(logcounts(sce))
# M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
# save(M, file="CTS_ShrinkM.RData", compress=TRUE) 
load("CTS_ShrinkM.RData")
dim(M) #96  96



############################################################
## 1) running BioTIP on downsized datasets, iteratively   ##
############################################################

######### setting parameters  ######################
localHVG.preselect.cut = 0.8 # A positive numeric value setting how many propotion of global HVG to be selected per cell cluster
getNetwork.cut.fdr = 0.2  # to construct RW network and extract co-expressed gene moduels
getTopMCI.gene.minsize = 10  # min number of genes in an identified CTS, 
                            # depending on the size of global HVG of the number of tested genes (e.g., here is 96 from scRT-PCR)
getTopMCI.n.states = 3  # A number setting the number of states to check the significance of DNB score (i.e. MCI) 
                        # This parameter can be enlarge when inputting more clusters
MCIbottom = 2:4  # A number setting the threshold to prioritize the initially selected CTS candidates.
               # In our experiment, a number between 2 and 4 generated expected resutls.

############### run BioTIP  ##################
n.scability <- 20
downsize <- 0.95  
n1 <- nrow(sce)
n2 <- ncol(sce)

for(i in 1:length(MCIbottom)){
  set.seed(102020)
  ## get n.scability sets of predictions
  res <- list()  
  flag = counts = 0
  repeat{ 
    cat(counts, '\t')
    counts = counts+1 # control the maxmum trials for searching two sets of signficant CTSs, given the current parameter settings
    sce.dnsize <- sce[sample(n1, floor(n1*downsize)), sample(n2, floor(n2*downsize))]
    samplesL <- split(rownames(colData(sce.dnsize)), f = sce.dnsize$Consensus.Cluster)
    this.run <- BioTIP.wrap(sce.dnsize, samplesL, subDir='stability', 
                            getTopMCI.n.states=getTopMCI.n.states, 
                            getTopMCI.gene.minsize= floor(getTopMCI.gene.minsize * downsize), 
                            MCIbottom=MCIbottom[i],
                            localHVG.preselect.cut=localHVG.preselect.cut, 
                            getNetwork.cut.fdr=getNetwork.cut.fdr, 
                            M=M,
                            permutation.method='both',
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
  save(res, flag, counts, file=paste0('stability/BioTIP_res_20runs_MCIbottom',MCIbottom[i],'.RData'))
  
}

###########################################################
## evaluate the stability of CT detection
## using the robustly detected 2 clusters (C9 and C10) as reference
###########################################################
df.CT =NULL
getTopMCI.gene.minsize
#[1] 10
for(m in 1:length(getTopMCI.gene.minsize)){
  CT.list <- list()
  Freq <- matrix(0, nrow=nlevels(sce$Consensus.Cluster), ncol=length(MCIbottom))#!!!
  rownames(Freq) <- levels(sce$Consensus.Cluster) #!!!
  colnames(Freq) <- MCIbottom 
  for(i in 1:length(MCIbottom)){
    load(file=paste0('stability/BioTIP_res_20runs_MCIbottom',MCIbottom[i],'.RData'))
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
  #      Frequency clusterID MCIbottom moduleSize
  # 1         0         1         2         10
  # 2        10        10         2         10
  # 3         0        11         2         10
  # 4         0        12         2         10
  
  df$clusterID <- factor(df$clusterID, levels =levels(sce$Consensus.Cluster)) #!!!
  df$moduleSize <- rep(getTopMCI.gene.minsize[m], nrow(df))
  
  df.CT <- rbind(df.CT, df)  
  save(CT.list, file= paste0('stability/CT.list_variable.Modulesize',getTopMCI.gene.minsize[m],'.Rdata'))
}
df.CT$Frequency = df.CT$Frequency/n.scability

save(df.CT, file= 'stability/CT.detection_variable.MCIbottom.Rdata')



##########################################################
## plot for variable MCIbottom, when modulesize==10     ##
##########################################################

pdf(file='stability/CT.detection_variable.MCIbottom.pdf')
ggplot(data=df, aes(x=clusterID, y=Frequency, group=MCIbottom)) +
  geom_line(aes(color=MCIbottom), size=1.5)+
  geom_point(aes(color=MCIbottom))
dev.off()

 
###########################################################
## evaluate the stability of CTS identification from downsampled datasets (n=20)
## for the C9 and C10 respectively
## by comparing to the original identification from the whole dataset (same parameters)
## and using size-controled random genes as a negative control
###########################################################
load('original_run/CTS.RData')
names(CTS)
#[1] "9"  "10" "9" 
C9.ref = CTS[which(names(CTS)=='9')]
C10.ref = CTS['10']

#load(file= 'stability/CT.detection_variable.MCIbottom.Rdata')
load(file='stability/CT.list_variable.Modulesize10.Rdata')
names(CT.list)
#[1] "2" "3" "4"

F1.C9 <- F1.C9.ctl <- F1.C10 <- F1.C10.ctl <- list()
for(i in 1:length(MCIbottom)){
  load(file=paste0('stability/BioTIP_res_20runs_MCIbottom',MCIbottom[i],'.RData'))
  F1.CTS.C9 <- F1.CTS.C10 <- array(dim=n.scability)
  F1.bk.C9 <- F1.bk.C10 <- array(dim=n.scability)
  for(j in 1:n.scability){
    # evaluate the F1 score for two CTS identificitons
    set2 <- res[[j]]$CTS.candidate[which(res[[j]]$significant)] 
    F1.CTS.C9[j] <- F_score.CTS(C9.ref, set2, weight=TRUE)$F1
    F1.CTS.C10[j] <- F_score.CTS(C10.ref, set2, weight=TRUE)$F1
   # negative control
    n.module <- lengths(set2)
    random.set2 <- lapply(n.module, function(x) sample(rownames(sce), x))
    F1.bk.C9[j] <- F_score.CTS(C9.ref, random.set2, weight=TRUE)$F1
    F1.bk.C10[j] <- F_score.CTS(C10.ref, random.set2, weight=TRUE)$F1
  } 
  F1.C9[[i]] =  F1.CTS.C9
  F1.C10[[i]] =  F1.CTS.C10
  F1.C9.ctl[[i]] =  F1.bk.C9
  F1.C10.ctl[[i]] =  F1.bk.C10
 }  
names(F1.C9) <- names(F1.C10) <- names(F1.C9.ctl) <- names(F1.C10.ctl) <- MCIbottom

save(F1.C9, F1.C10, F1.C9.ctl, F1.C10.ctl,
     file='stability/F1.CTS_MCIbottom.RData')

lapply(F1.C9, mean) %>% unlist()
#[1]0.3328190 0.3252816 0.4289059
lapply(F1.C9.ctl, mean) %>% unlist()
# 0.08125084 0.08401564 0.08445800 
lapply(F1.C10, mean) %>% unlist()
#0.3723308 0.3258299 0.2002311 
lapply(F1.C10.ctl, mean) %>% unlist()
#0.04866308 0.04572070 0.02939498 

df <- data.frame(F1 = c(unlist(F1.C9.ctl),unlist(F1.C9),
                        unlist(F1.C10.ctl),unlist(F1.C10)),
                 MCIbottom = c(rep(MCIbottom, lengths(F1.C9.ctl)),
                               rep(MCIbottom, lengths(F1.C9)), 
                               rep(MCIbottom, lengths(F1.C10.ctl)),
                               rep(MCIbottom, lengths(F1.C10))),
                 CT = c(rep('C9', sum(lengths(F1.C9.ctl))+sum(lengths(F1.C9))), 
                        rep('C10', sum(lengths(F1.C10.ctl))+sum(lengths(F1.C10)))),
                CTS = c(rep('ctl', sum(lengths(F1.C9.ctl))), rep('CTS', sum(lengths(F1.C9))), 
                        rep('ctl', sum(lengths(F1.C10.ctl))), rep('CTS', sum(lengths(F1.C10)))) 
)
df$MCIbottom <- as.character(df$MCIbottom) 
df$CT <- factor(df$CT, levels = c('C9','C10')) 
df$CTS <- factor(df$CTS, levels = c('ctl','CTS')) 
df$color.lab <- paste0(df$MCIbottom,df$CTS)
df$color.lab <- factor(df$color.lab, levels=c('2ctl','2CTS','3ctl','3CTS','4ctl','4CTS'))

# To ensure that easy and difficult CT state have equal influence on the final score, 
# we normalized the scores at each CT state across the different MCIbottom setting.
df$Norm.F1 <- Normalize.F1(df$F1)
 
library(ggbeeswarm)
library(grid)
library(gridExtra)
library(ggpubr)

#my_comparisons <- list( c("2ctl", "2CTS"), c("3ctl", "3CTS"), c("4ctl", "4CTS") )

p1 <- ggplot(data=df, aes(x=CT, y=F1, color=color.lab)) +
  geom_violin(position = position_dodge(width = 0.5)) +
  geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
   xlab('predicted for cluster') + ylab('F1 score for the CTS') +
  scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue")) 

p2 <- ggplot(data=df, aes(x=CT, y=F1, color=color.lab)) +
  geom_violin(position = position_dodge(width = 0.5)) +
  geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
  xlab('predicted for cluster') + ylab('F1 score for the CTS') +
  scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue")) 

pdf(file='stability/CTS.F1_variable.MCIbottom.pdf')
gridExtra::grid.arrange(
   p1, p2, top = "CTS stability", bottom = "different MCIbottom settings"
  ,ncol=1)
dev.off()

tmp <- compare_means(Norm.F1 ~ color.lab, data=df, group='CT', method='t.test')
subset(tmp, group1=='2ctl' & group2=='2CTS')
# # A tibble: 2 x 9
# CT    .y.   group1 group2       p p.adj p.format p.signif method
# <fct> <chr> <chr>  <chr>    <dbl> <dbl> <chr>    <chr>    <chr> 
#1 C9    Norm.F1 2ctl   2CTS   0.00157 0.039 0.00157  **       T-test
#2 C10   Norm.F1 2ctl   2CTS   0.00201 0.039 0.00201  **       T-test
subset(tmp, group1=='3ctl' & group2=='3CTS')
# # A tibble: 2 x 9
# CT    .y.   group1 group2       p p.adj p.format p.signif method
# <fct> <chr> <chr>  <chr>    <dbl> <dbl> <chr>    <chr>    <chr> 
#1 C9    Norm.F1 3ctl   3CTS   0.00165 0.039 0.00165  **       T-test
#2 C10   Norm.F1 3ctl   3CTS   0.0111  0.19  0.01114  *        T-test
subset(tmp, group1=='4ctl' & group2=='4CTS')
# # A tibble: 2 x 9
# CT    .y.   group1 group2        p  p.adj p.format p.signif method
# <fct> <chr> <chr>  <chr>     <dbl>  <dbl> <chr>    <chr>    <chr> 
# 1 C9    Norm.F1 4ctl   4CTS   0.000271 0.0079 0.00027  ***      T-test
# 2 C10   Norm.F1 4ctl   4CTS   0.0627   0.94   0.06273  ns       T-test

 