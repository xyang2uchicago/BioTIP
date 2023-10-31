# This code doing:
# 1) Annotate cell identities by comparing to the original cluster IDs.
#
# 2) Summarize the CT identifications
#
# 3) Calculating Jaccard-similarity for CT detections
#
# 4) Calculate stability score between the original CTS identification (Author-published consensus clusters) 
## to the CTS identification of each set of cell cluster identities.
#
# 5) Verify if the CTS identified KIT, the experimetnally validated early-warning singals for the PS->CM transition
#
# 6) From the predictions have identified KIT, we generate a proxy 'golden standard' of CTS genes,
## which were identifed by at least 80% of the selected methods.
## Then calculate Prcision Recall and AUC for each prediction, and plot ROC plot
# 
# last update 1/27/2022
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
library(viridis)
library(gridExtra)  # required by grid.arrange()
library(ggpubr)  # required by p-values on bar
library(ROCR)


####################################################
##                                                ##
## load the R object                              ##                         
## focusing on 1.3k cells of the HEP processes    ##
## refer to BioTIP_E8.25_mesoderm_add.clusters.R  ##
##                                                ##
####################################################

load('sce_hESC.RData')
sce
# class: SingleCellExperiment 
# dim: 96 929

# logmat <- as.matrix(logcounts(sce))
# M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
# save(M, file="CTS_ShrinkM.RData", compress=TRUE) 
# load(file="CTS_ShrinkM.RData")
# dim(M) #96  96



######################################################################
##                                                                  ##
## 1) Summarize the CT identifications for                          ##
## the three of the original predicted CT states which are          ## 
## later HEP C15; Endothelium C6, and early HEP S13                 ##
##                                                                  ##
######################################################################

running.methods <- list.dirs()
running.methods <- running.methods[-1]
methods <- lapply(running.methods, function(x) unlist(strsplit(x, "/C_"))[2]) %>% unlist()
x <- which(is.na(methods))
if(length(x)>0){
  running.methods = running.methods[-x]
  methods = methods[-x]
}

CT <- list() 
for(i in 1:length(running.methods)){
  load(file=paste0(running.methods[i],'/BioTIP.res.RData'))
  if(length(res)>1) {
    CT[[i]] <- res$significant[which(res$significant)] %>% names() %>% unique()
  } else CT[[i]] <- NA
}
names(CT) <- paste0("C_", methods)

## manually add teh QuanTC's prediction about CT states,
## See Figure S2 for details
CT$QuanTC_run <- c('C3','C2','C4')
CT

###########################################################################
##                                                                       ##
## 1) Annotate cell identities by comparing to the original cluster IDs. ##
##                                                                       ##
###########################################################################
colnames(colData(sce))
# [1] "CollectionTime"    "Consensus.Cluster" "Phenotype.FACS"   
# [4] "mesodermPath"      "C_SNNGraph_k10"    "C_SNNGraph_k20"   
# [7] "C_Soft"            "C_consensus_ks4"   "C_consensus_ks5"  
# [10] "C_consensus_ks6"   "C_consensus_ks7"   "C_consensus_ks8"  
# [13] "C_consensus_ks9"   "C_consensus_ks10"  "C_Leiden_0.4"     
# [16] "C_Leiden_0.8"      "C_Leiden_1.2" 
x <- grep("C_", colnames(colData(sce)))
colnames(colData(sce))[x]
 
best.matched.clusterID <- data.frame()
query = c('9','10','7','1') # the cluster ID for bifurcation at PS, CM, EP, and Ecoderm
for(i in x){
  tmp <- table(sce$Consensus.Cluster, colData(sce)[,i])
  n1 <- nrow(tmp)
  n2 <- ncol(tmp)
  Jaccard_cell = tmp*0
  for(a in 1:n1){
    for(b in 1:n2){
      Jaccard_cell[a,b] <- tmp[a,b]/(sum(tmp[a,])+sum(tmp[,b])-tmp[a,b])
    }
  }
  best.match <- colnames(Jaccard_cell)[apply(Jaccard_cell, 1, which.max)]  
  names(best.match) <- rownames(Jaccard_cell)
  best.matched.clusterID <- rbind(best.matched.clusterID, 
                                  best.match[query])
}
rownames(best.matched.clusterID) <- colnames(colData(sce))[x]
colnames(best.matched.clusterID) <- c('PS','CP','EP','EctP')

## additionally recalculate without the TC among QuanTC's clusters
i= grep('Soft', colnames(colData(sce))) # 7
tmp <- table(sce$Consensus.Cluster, colData(sce)[,i])
tmp = tmp[,1:(ncol(tmp)-1)] # masking TC 
n1 <- nrow(tmp)
n2 <- ncol(tmp)
F1 = tmp*0
for(a in 1:n1){
  for(b in 1:n2){
    F1[a,b] <- tmp[a,b]/(sum(tmp[a,])+sum(tmp[,b])-tmp[a,b])
  }
}
best.match <- colnames(F1)[apply(F1, 1, which.max)]  
names(best.match) <- rownames(F1)
best.matched.clusterID <- rbind(best.matched.clusterID, 
                                'C_Soft.wo.TC'= best.match[query])

### add the best matches for the additional BioTIP runs and the original QuanTC run
setdiff(names(CT), rownames(best.matched.clusterID))
# [1] "C_CollectionTime"          "C_CollectionTime_allcells"
# [3] "C_Consensus_929cells"      "C_Consensus_allcells"     
# [5] "QuanTC_run"

# Now, manually add the best-match infor for these additional BioTIP runs
n <- nrow(best.matched.clusterID) # 8
best.matched.clusterID <- rbind(best.matched.clusterID, 
                                c('Day2.5', 'Day3', 'Day2.5','Day1'),c('Day2.5', 'Day3', 'Day2.5','Day1'), 
                                query, query,
                                c('C3','C2','C4','C1')) 
rownames(best.matched.clusterID)[(n+1):(n+5)] <- setdiff(names(CT), rownames(best.matched.clusterID))


save(best.matched.clusterID, file='best.matched.clusterID.RData')
best.matched.clusterID

################################################################################################
##                                                                                            ##
## 3) calculate Jaccard-similarity score for CT identifications,                              ##
## between each predictions and that of BioTIP using the published consensus clusters         ##
## after finding the reference cluster ID (i.e. best-matched IDs) for the QuanTC clusters     ##
## And maually adding the 'best matches' for additional BioTIP/QuanTC runs of this datasets   ##
##                                                                                            ##
################################################################################################

jaccard.CT <- array()
for(i in 1:length(CT)){
  x <- names(CT)[i]
  # in case BioTIP analysis of 929 cells excluding endothelia progenitors, goldstandard is the first two reference 
  n <- ifelse(grepl('allcells', x) | grepl('QuanTC_run', x), 3, 2)  
  reference <- best.matched.clusterID[x,1:n]
  jaccard.CT[i] <- jaccard.sim(reference, CT[[x]]) 
}
names(jaccard.CT) <- names(CT)
round(jaccard.CT, 2)


################################################################
##                                                            ##
## 4) evaluate the F1 score for the identified CTS genes      ##
## for any commonly identifed CTs from two identifications    ##
##                                                            ##
################################################################
# load the orignial identification of three CTSs along the PS->CM trajectory
# refer to data_matrix_Cluster_CMPath.R 
# load('../Cluster/CTS_maxMCI.RData')
# CTS <- CTS[c(1,2,8)] # this line is particular for the old saved R object which had all candidates saved
# save(CTS, file='original_run/CTS.RData')
load('original_run/CTS.RData')
lengths(CTS) 
#  9 10  9 
# 18 19 23 
CTS.reference <- CTS
rm(CTS)

# load the CTS prediction with different clustering methods
running.methods <- list.dirs()
running.methods <- running.methods[-1]
methods <- lapply(running.methods, function(x) unlist(strsplit(x, "./C_"))[2]) %>% unlist()
x <- which(is.na(methods))
if(length(x)>0) { 
  methods <- methods[-x]
  running.methods <- running.methods[-x]
  }
## manually adding the QuanTC's predictions
running.methods <- c(running.methods, './QuanTC_run')
methods <- c(methods,'QuanTC_run')

# review QuanTC's identification of transition genes
load('QuanTC_run/CTS.RData')
CTS
# $C3.to.C2
# [1] "GATA6"  "BMP2"   "DKK1"   "EOMES"  "GSC"    "PDGFRA" "GATA4"  "EVX1"  
# [9] "WNT5A"  "FOXC1"  "FZD7"   "MIXL1" 
# 
# $C2.to.C4
# [1] "SOX17"

res.PS <- get.multiple.F1.score(match.col.name=c('PS'),
                                CTS.reference = unique(CTS.reference[[1]],CTS.reference[[3]]),
                                QuanTC.genes.list.id = 1, ## using the published 9 gene (Fig for QuanTC results
                                running.methods=running.methods, 
                                methods=methods,
                                sce=sce)
res.CP <- get.multiple.F1.score(match.col.name=c('CP'),
                                CTS.reference = CTS.reference[['10']],
                                QuanTC.genes.list.id = 2,
                                running.methods=running.methods, 
                                methods=methods,
                                sce=sce)

x <- which(res.PS$F1.scores$PS>0) 
t.test(res.PS$F1.scores$PS[x], res.PS$F1.ctl$PS[x])
# t = 1.9847, df = 10.862, p-value = 0.07301
x <- which(res.CP$F1.scores$CP>0) 
t.test(res.CP$F1.scores$CP[x], res.CP$F1.ctl$CP[x])
#t = 5.1876, df = 4.1623, p-value = 0.005892

F1.PS <- res.PS$F1.scores$PS 
F1.CP <- res.CP$F1.scores$CP 
F1.PS.ctl <- res.PS$F1.ctl$PS 
F1.CP.ctl <- res.CP$F1.ctl$CP 

## Regarding F1 scores for the gene module with the highest DNB score
# which is C9_19g of the original run
F1.DNB.PS <- get.multiple.F1.score(match.col.name=c('PS'),
                                   CTS.reference = CTS.reference[[1]],
                                   QuanTC.genes.list.id = 1, ## using the published 9 gene (Fig for QuanTC results
                                   running.methods=running.methods, 
                                   methods=methods,
                                   sce=sce,
                                   best.matched.clusterID=best.matched.clusterID)
t.test(F1.DNB.PS$F1.scores$PS, F1.DNB.PS$F1.ctl$PS)  # p-value = 0.1277 #!!!!!
F1.PS.DNB <- F1.DNB.PS$F1.scores$PS
F1.PS.ctl.DNB <- F1.DNB.PS$F1.ctl$PS


nn = length(F1.PS.ctl)  # 19
Normalized.F1 <- Normalize.F1(c(F1.PS.ctl, F1.PS, F1.CP.ctl, F1.CP, F1.PS.ctl.DNB, F1.PS.DNB))
Norm.F1.PS.ctl <- Normalized.F1[1:nn] 
Norm.F1.PS <- Normalized.F1[(nn+1):(2*nn)] 
Norm.F1.CP.ctl <- Normalized.F1[(2*nn+1):(3*nn)]
Norm.F1.CP <- Normalized.F1[(3*nn+1):(4*nn)] 

Norm.F1.PS.ctl.DNB <- Normalized.F1[(4*nn+1):(5*nn)]
Norm.F1.PS.DNB <- Normalized.F1[(5*nn+1):(6*nn)] 

###################################################################################################################
##                                                                                                               ##
# 5) Verify if the CTS identified the estabilished early-warning singals which are                               ##
## The validated endoderm-mesoderm branching marker KIT (Bargaje2017 Fig3)                                       ##
##                                                                                                               ##
###################################################################################################################
# load the CTS prediction for the PS state with different clustering methods
PS.marker <-  c('KIT','DKK1')  #c('KIT', 'DKK1', 'PDGFRA')   
  
PS.CTS <- NULL 
for(i in 1:length(running.methods)){
  if(methods[i] == 'QuanTC_run') {
    load(file=paste0(running.methods[i],'/CTS.RData'))
    set2 <- CTS 
    names(set2) <- c('C3','C4') # normalzie the name in order to use the weight of F1 score
    best.match <- best.matched.clusterID[methods[i],'PS']
    x <- grep(best.match, names(set2))
    if(length(x)>0) tmp <- PS.marker %in% unlist(set2[x]) else tmp <- rep('FALSE', length(PS.marker))
    PS.CTS = rbind(PS.CTS, tmp)
    } else {
      load(file=paste0(running.methods[i],'/BioTIP.res.RData'))
      if(length(res)>1){ # flag 1
        names(res)
        #[1] "CTS.candidate" "CTS.score"     "Ic.shrink"     "significant" 
        set2 <- res$CTS.candidate[which(res$significant)]
      
        best.match <- best.matched.clusterID[paste0('C_',methods[i]),'PS']
        x <- grep(best.match, names(set2))
        if(length(x)>0) tmp <- PS.marker %in% unlist(set2[x]) else tmp <- rep('FALSE', length(PS.marker))
        PS.CTS = rbind(PS.CTS, tmp)
       } else  { # flag 1
        PS.CTS <- rbind(PS.CTS, rep('FALSE', length(PS.marker)))
      }
    }
}
rownames(PS.CTS) <- methods
colnames(PS.CTS) <- PS.marker

save(PS.CTS, file='PS.CTS.RData')




############################################################################################
##                                                                                         
## # 6) Additionally, from the predictions have identified KIT successfully for the PS->CM transition,   
## we generate a proxy 'golden standard' of CTS genes,                                     
## which were identifed by at least 4 out of 5 of the selected methods.                           
##                                                                                         
############################################################################################
## select the clustering methods with an optimal parameter settings that meet the following creteial
## a) must detect both KIT or DKK1
x <- apply(PS.CTS, 1, function(x) ifelse(x,1,0) %>% sum())
(select <- names(x)[which(x>=2)])
# [1] "Consensus_929cells" "consensus_ks10"     "Leiden_0.4"         "SNNGraph_k10"       "Soft"    
(n.select <- floor(length(select) * 0.8)) # 4

GS.CTS = NULL
for(i in 1:length(select)){
  if(select[i] == 'QuanTC_run') {
    load(file=paste0(select[i],'/CTS.RData'))
    set2 <- CTS[[1]]
    GS.CTS = c(GS.CTS, unlist(set2[x]))
     } else  {
    load(file=paste0('C_',select[i],'/BioTIP.res.RData'))
    set2 <- res$CTS.candidate[which(res$significant)]
    best.match <- best.matched.clusterID[paste0('C_',select[i]),'PS']
    x <- grep(best.match, names(set2))
    GS.CTS = c(GS.CTS, unlist(set2[x]))
     }
}
x <- table(GS.CTS)
GS.CTS <- names(x)[which(x>= n.select)]
GS.CTS
#[1] "DKK1"  "DLL3"  "EVX1"  "GATA6" "GSC"   "KDR"   "KIT"   "MIXL1" "T"   
(n.GS <- length(GS.CTS))  # 9
write.table(GS.CTS, file=paste0('GS.CTS.PS_at.least.',n.select,'identificaitons_fr.',length(select),'predicitons.txt'),
            row.names=F, col.names=F)

## generate dataframe for ROC plotting and calcuate the Precision, Recall
df <- data.frame(gene=rownames(sce), 
                 GS.CTS= ifelse(rownames(sce) %in% GS.CTS, 1,0))

PS.Recall <- PS.Precision <- NULL 
for(i in 1:length(select)){
  if(select[i] == 'QuanTC_run') {
    load(file=paste0(select[i],'/CTS.RData'))
    set2 <- CTS 
    names(set2) <- c('C3','C4') # normalzie the name in order to use the weight of F1 score
    best.match <- best.matched.clusterID[select[i],'PS']
    } else {
      load(file=paste0('C_',select[i],'/BioTIP.res.RData'))
      set2 <- res$CTS.candidate[which(res$significant)]
      best.match <- best.matched.clusterID[paste0('C_',select[i]),'PS']
   } 
  x <- grep(best.match, names(set2))
  PS.Precision <- c(PS.Precision, lapply(set2[x], function(y) length(which(GS.CTS %in% y))/length(unlist(y))) %>% 
                      unlist() %>% 
                      mean())
  PS.Recall <- c(PS.Recall, lapply(set2[x], function(y) length(which(GS.CTS %in% y))/n.GS) %>% 
                      unlist() %>% 
                      mean()) 
  df <- cbind(df, ifelse(rownames(sce) %in% unlist(set2), 1, 0))
  }

names(PS.Recall) <- names(PS.Precision) <- select
colnames(df) <- c('Gene','PS.GS', select)

# get the auc value
auc <- NULL
for(i in 1:length(select)){
  pred <- ROCR::prediction(df[,(2+i)], df$PS.GS)
  auc_ROCR <- performance(pred, measure = "auc")
  auc[i] <- auc_ROCR@y.values[[1]]
}  
names(auc) <- select
round(auc, 2)
# Consensus_929cells     consensus_ks10         Leiden_0.4 
# 0.89               0.94               0.94 
# SNNGraph_k10               Soft 
# 0.69               0.74 

##########################################
## reproting table #######################
tb <- data.frame(Jaccard.CT = round(jaccard.CT,3),
                 F1.PS.ctl = round(F1.PS.ctl,3),
                 F1.PS = round(F1.PS,3),
                 F1.CP.ctl = round(F1.CP.ctl,3),
                 F1.CP = round(F1.CP,3),
                 Norm.F1.PS.ctl = round(Norm.F1.PS.ctl,3),
                 Norm.F1.PS = round(Norm.F1.PS,3),
                 Norm.F1.CP.ctl = round(Norm.F1.CP.ctl,3),
                 Norm.F1.CP = round(Norm.F1.CP,3),
                 F1.PS.ctl.DNB = round(F1.PS.ctl.DNB,3),
                 F1.PS.DNB = round(F1.PS.DNB,3),
                 Norm.F1.PS.ctl.DNB = round(Norm.F1.PS.ctl.DNB,3),
                 Norm.F1.PS.DNB = round(Norm.F1.PS.DNB,3)
)
tb <- cbind(tb, PS.CTS)
tb$Recall <- tb$Precision <- tb$AUC <- rep(NA, nrow(tb))
tb[paste0('C_',select),'Recall'] <- PS.Recall
tb[paste0('C_',select),'Precision'] <- PS.Precision
tb[paste0('C_',select),'AUC'] <- auc
tb$detected.CT <- lapply(CT, toString) %>% unlist()
tb$best.matched.cluster <- apply(best.matched.clusterID,1, toString)
nrow(tb) #[1] 19

# reorder to show
myorder <- c('C_Consensus_929cells', 'C_Consensus_allcells',
             'C_CollectionTime','C_CollectionTime_allcells',
             paste0('C_consensus_ks',4:10),
             'C_SNNGraph_k10','C_SNNGraph_k20', 
             paste0('C_Leiden_',c(0.4, 0.8, 1.2)),
             'C_Soft.wo.TC', 'C_Soft','QuanTC_run' 
             )
  
(mytable <- tb[myorder,
   c('best.matched.cluster','detected.CT','Jaccard.CT', 
     'F1.PS.ctl', 'Norm.F1.PS.ctl', 'F1.PS', 'Norm.F1.PS',  
     'F1.CP.ctl', 'Norm.F1.CP.ctl', 'F1.CP', 'Norm.F1.CP',
     'KIT',  'DKK1', 'Precision',  'Recall', 'AUC', 
     'F1.PS.ctl.DNB', 'Norm.F1.PS.ctl.DNB', 'F1.PS.DNB','Norm.F1.PS.DNB')])
write.table(mytable,   file='jaccard.CT_F1.CTS.txt', 
            sep='\t')


###############
## PLOT #######
###############
## generate ROC plot ##########################
for(i in 1:length(select)){
  pred <- ROCR::prediction(df[,(2+i)], df$PS.GS)
  perf <- ROCR::performance(pred,"tpr","fpr")
  add.line = ifelse(i==1, FALSE, TRUE)
  ROCR::plot(perf,colorize=FALSE, add=add.line, col=i,
             main='hESC Bargaje2017 data')
}
legend(0.5,0.5, select, text.col=1:length(select), bty='n')
dev.copy2pdf(file='ROC.PS.CTS_clustering.inputs_BioTIP.outputs.pdf')


## generate bar plot #######################
df <- read.table('jaccard.CT_F1.CTS.txt', header=T, sep='\t')
tmp <- rownames(df)
tmp = gsub('C_', 'BioTIP_', tmp)
df$method = gsub('QuanTBioTIP_run', 'QuanTC_Soft', tmp)

## select certain rows to plot 
df <- df[c("QuanTC_run", "C_CollectionTime", "C_Soft.wo.TC", #"C_CollectionTime_allcells", "C_Consensus_allcells" ,
           "C_Soft", 
           "C_Consensus_929cells", "C_consensus_ks7","C_consensus_ks8","C_consensus_ks10",
           "C_SNNGraph_k10", "C_SNNGraph_k20",
           "C_Leiden_0.4","C_Leiden_0.8" ,"C_Leiden_1.2", ),]
                    
df$method=factor(df$method, levels=df$method)

p1 <- ggbarplot(df, "method", "Jaccard.CT", orientation = "horiz") +
     geom_col(aes(fill = Jaccard.CT)) + 
     scale_fill_gradient2(low = "white",high = "green") +
     theme_bw(base_size=10)  +        # use a black-and-white theme with set font size
    theme(legend.position = "none") + coord_flip()   
p2 <- ggbarplot(df, "method", "Norm.F1.PS", orientation = "horiz") +
     geom_col(aes(fill = Norm.F1.PS)) + 
     scale_fill_gradient2(low = "white",high = "blue") +
     theme_bw(base_size=10)    +      # use a black-and-white theme with set font size
     theme(legend.position = "none")  + coord_flip()  
p3 <- ggbarplot(df, "method", "Norm.F1.CP", orientation = "horiz") +
  geom_col(aes(fill = Norm.F1.CP)) + 
  scale_fill_gradient2(low = "white",high = "blue") +
  theme_bw(base_size=10)    +      # use a black-and-white theme with set font size
  theme(legend.position = "none")  + coord_flip()  
p4 <- ggbarplot(df, "method", "KIT", orientation = "horiz") +
     geom_col(aes(fill=KIT)) + 
     scale_fill_manual(values=c("grey", "dark red")) +
     theme_bw(base_size=10) + 
     theme(legend.position = "none")   + coord_flip()      
p5 <- ggbarplot(df, "method", "DKK1", orientation = "horiz") +
  geom_col(aes(fill=DKK1)) + 
  scale_fill_manual(values=c("grey", "dark red")) +
  theme_bw(base_size=10) + 
  theme(legend.position = "none")   + coord_flip()      

pdf(file='bar_robustness.pdf', width=15, height=3)
gridExtra::grid.arrange(
     p1, p2, p3, p4,p5 
     ,ncol=5)
dev.off()

