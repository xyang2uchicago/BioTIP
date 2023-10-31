# This code doing:
# 1) Run BioTIP given different outputs of different clustering methods.
#
# last update 1/26/2022
# by Holly yang

setwd('F:/projects/BioTIP/result/Bargaje2017_EOMES/hESC_Bargaje2017_robustness/')

# the following funciton has been updated on github, but before updating the package to Bioconductor, you need to call it locally to test
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/ForJennifer/optimize.sd_selection.R') 
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

logmat <- as.matrix(logcounts(sce))
M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
save(M, file="CTS_ShrinkM.RData", compress=TRUE) 
dim(M) #96  96


############################################################
## 1) running BioTIP on different cell clustering outputs ##
############################################################

######### setting parameters  ######################
localHVG.preselect.cut = 0.8 # A positive numeric value setting how many propotion of global HVG to be selected per cell cluster
getNetwork.cut.fdr = 0.2  # to construct RW network and extract co-expressed gene moduels
getTopMCI.gene.minsize = 10  # min number of genes in a identified CTS
getTopMCI.n.states = 3  # A number setting the number of states to check the significance of DNB score (i.e. MCI) 
                         # This parameter can be enlarge when inputting more clusters, usually the expected number of transition states +1
shuffle.both = TRUE
MCIbottom = 3  # A number setting the threshold to prioritize the initially selected CTS candidates.
               # In our experiment, a number between 2 and 4 generated expected resutls.

############### all new clustering outputs  ##################
x <- grep("C_", colnames(colData(sce)))
colnames(colData(sce))[x]
# [1] "C_SNNGraph_k10"   "C_SNNGraph_k20"   "C_Soft"           "C_consensus_ks4" 
# [5] "C_consensus_ks5"  "C_consensus_ks6"  "C_consensus_ks7"  "C_consensus_ks8" 
# [9] "C_consensus_ks9"  "C_consensus_ks10" "C_Leiden_0.4"     "C_Leiden_0.8"    
# [13] "C_Leiden_1.2"  

# set.seed(102020)
for(i in 1:length(x)){
  samplesL <- split(rownames(colData(sce)), f = colData(sce)[x[i]])
  (tmp = lengths(samplesL))

#  getTopMCI.n.states = ifelse(length(samplesL)<10, 3, 4)
  
  res <- BioTIP.wrap(sce, samplesL, subDir=colnames(colData(sce))[x[i]],
                     getTopMCI.n.states=getTopMCI.n.states, 
                     getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                     MCIbottom=MCIbottom,
                     localHVG.preselect.cut=localHVG.preselect.cut, 
                     getNetwork.cut.fdr=getNetwork.cut.fdr, 
                     M=M,
                     permutation.method='both',
                     verbose=TRUE, plot=TRUE)           
  save(res, file=paste0(colnames(colData(sce))[x[i]],'/BioTIP.res.RData'))  
  
}

# # If you want to keep the simulation resutls but tune the curoff of signficance, 
# # e.g. we tried empirical.IC.p.cut=0.05 and  local.IC.p=TRUE here
# source("F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/BioTIP.tune.parameter.R")
# for(i in 11:length(x)){
#   samplesL <- split(rownames(colData(sce)), f = colData(sce)[x[i]])
#   (tmp = lengths(samplesL))
# 
#   #  getTopMCI.n.states = ifelse(length(samplesL)<10, 3, 4)
# 
#   res <- BioTIP.tune.parameter(sce, samplesL, subDir=colnames(colData(sce))[x[i]],globa.HVG.select=FALSE,
#                                # dec.pois=dec.pois, n.getTopHVGs=4000,
#                                getTopMCI.n.states=getTopMCI.n.states, getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
#                                n.getMaxMCImember = 2, MCIbottom=MCIbottom,
#                                local.HVG.optimize =TRUE, localHVG.preselect.cut=localHVG.preselect.cut, localHVG.runs=100,
#                                getNetwork.cut.fdr=getNetwork.cut.fdr, 
#                                n.permutation = 100,  M=M,
#                                permutation.method='both',
#                                plot=TRUE) 
#   #save(res, file=paste0(colnames(colData(sce))[x[i]],'/BioTIP.res.RData'))
# }

## Additional run for soft thresholding with TC
samplesL <- split(rownames(colData(sce)), f = colData(sce)$C_Soft)
(tmp = lengths(samplesL))
# C1  C2  C3  C4  TC 
# 231 153 249 104 192
samplesL = samplesL[1:(length(samplesL)-1)]

res <- BioTIP.wrap(sce, samplesL, subDir='C_Soft.wo.TC',
                   getTopMCI.n.states=getTopMCI.n.states, 
                   getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                   MCIbottom=MCIbottom,
                   localHVG.preselect.cut=localHVG.preselect.cut, 
                   getNetwork.cut.fdr=getNetwork.cut.fdr, 
                   M=M,
                   permutation.method='both',
                   verbose=TRUE, plot=TRUE)           
save(res, file='C_Soft.wo.TC/BioTIP.res.RData')  

## additional run for time-course day0-3
samplesL <- split(rownames(colData(sce)), f = colData(sce)$CollectionTime)
(tmp = lengths(samplesL))
# Day0   Day1 Day1.5   Day2 Day2.5   Day3   Day4   Day5 
# 231    166     93    211    122    106      0      0

res <- BioTIP.wrap(sce, samplesL, subDir='C_CollectionTime',
                   getTopMCI.n.states=getTopMCI.n.states, 
                   getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                   MCIbottom=MCIbottom,
                   localHVG.preselect.cut=localHVG.preselect.cut, 
                   getNetwork.cut.fdr=getNetwork.cut.fdr, 
                   M=M,
                   permutation.method='both',
                   verbose=TRUE, plot=TRUE)           
save(res, file='C_CollectionTime/BioTIP.res.RData')  


## additional run for all collected time-course cells, day0-5 
load('F:/projects/BioTIP/result/Bargaje2017_EOMES/sce.RData')
sce
# class: SingleCellExperiment 
# dim: 96 1896 
samplesL <- split(rownames(colData(sce)), f = colData(sce)$CollectionTime)
(tmp = lengths(samplesL))
# Day0   Day1 Day1.5   Day2 Day2.5   Day3   Day4   Day5 
#  231    166     93    211    122    258    456    359 

res <- BioTIP.wrap(sce, samplesL, subDir='C_CollectionTime_allcells',
                   getTopMCI.n.states=getTopMCI.n.states, 
                   getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                   MCIbottom=MCIbottom,
                   localHVG.preselect.cut=localHVG.preselect.cut, 
                   getNetwork.cut.fdr=getNetwork.cut.fdr, 
                   M=M,
                   permutation.method='both',
                   verbose=TRUE, plot=TRUE)           
save(res, file='C_CollectionTime_allcells/BioTIP.res.RData')  

## additional run for all published consensus clusters
samplesL <- split(rownames(colData(sce)), f = colData(sce)$Consensus.Cluster)
(tmp = lengths(samplesL))
#   1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
#  77 107 130 299 157  64  70  81  67  77 161  86  86  80  37 173  95  49

res <- BioTIP.wrap(sce, samplesL, subDir='C_Consensus_allcells',
                   getTopMCI.n.states=getTopMCI.n.states, 
                   getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                   MCIbottom=MCIbottom,
                   localHVG.preselect.cut=localHVG.preselect.cut, 
                   getNetwork.cut.fdr=getNetwork.cut.fdr, 
                   M=M,
                   permutation.method='both',
                   verbose=TRUE, plot=TRUE)           
save(res, file='C_Consensus_allcells/BioTIP.res.RData')  

## additionally re-run for the published consensus clusters, 929 cells, 
## compared to the publication, we decrease the initial cutoff of CTS score (in the DNB model) from 4 to 2 to scan ore CTS candidates
load('sce_hESC.RData')
sce
# class: SingleCellExperiment 
# dim: 96 929
samplesL <- split(rownames(colData(sce)), f = colData(sce)$Consensus.Cluster)
(tmp = lengths(samplesL))

res <- BioTIP.wrap(sce, samplesL, subDir='C_Consensus_929cells',
                   getTopMCI.n.states=getTopMCI.n.states, 
                   getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                   MCIbottom=MCIbottom,
                   localHVG.preselect.cut=localHVG.preselect.cut, 
                   getNetwork.cut.fdr=getNetwork.cut.fdr, 
                   M=M,
                   permutation.method='both',
                   verbose=TRUE, plot=TRUE)             
save(res, file='C_Consensus_929cells/BioTIP.res.RData')  


## copy QuanTC-identified transition genes
# BEGIN DO not repeat read and save QuanTC-predicted CTS genes   ##
CTS <- list()
tmp <- read.table('F:/projects/BioTIP/result/Bargaje2017_EOMES/QuanTC-modified/Output/k4_4_time/transition_gene_name1.1.txt',
                  header=FALSE)
CTS[['C3.to.C2']] <- tmp[,1]
tmp <- read.table('F:/projects/BioTIP/result/Bargaje2017_EOMES/QuanTC-modified/Output/k4_4_time/transition_gene_name1.2.txt',
                  header=FALSE)
CTS[['C3.to.C2']] <- c(QuanTC[['C3.to.C2']], tmp[,1])

tmp <- read.table('F:/projects/BioTIP/result/Bargaje2017_EOMES/QuanTC-modified/Output/k4_4_time/transition_gene_name2.1.txt',
                  header=FALSE)
CTS[['C2.to.C4']] <- tmp[,1]

CTS
# [['C3.to.C2]]
# [1] "GATA6"  "BMP2"   "DKK1"   "EOMES"  "GSC"    "PDGFRA" "GATA4"  "EVX1"
# [9] "WNT5A"  "FOXC1"  "FZD7"   "MIXL1" "LEFTY1"
#
# [['C2.to.C4']]
# [1] "SOX17"
save(CTS, file='QuanTC_run/CTS.RData')
# END DO not repeat read and save QuanTC-predicted CTS genes   ##
###################################################################

