# This code doing:
# 1) Annotate cell identities by comparing to the original cluster IDs.
#
# 2) Summarize the CT identifications
#
# 3) Calculating Jaccard-similarity for CT detections
## Note that three of the oribinal predicted CT states (later HEP C15; Endothelium C6, and early HEP S13) 
## are considered when teting on the 1362 cells cover the major HEP processes 
## (i.e., cells of the original C3, C13, C10, C6, C15, C7).
## One originally predicted CT state (Muscle/Mesenchyme C16) was excluded because it has been excluded from the 1362 tested cells. 
#
# 4) Calculate stability score between the original CTS identification (Author-published consensus clusters) 
## to the CTS identification of each set of cell cluster identities.
#
# 5) Verify if the CTS identified ETV2, the established marker for the HEP bifurcation
#
# 6) Due to no gold standard of CTS, we skip the step to estiamte the 
## Prcision Recall and AUC for each prediction, and plot ROC plot
# 
# last update 2/2/2022
# by Holly yang
setwd('F:/projects/scRNA/results/GSE87038_gastrulation/uncorrected/E8.25_mesoderm_robustness/')
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

load('sce_E8.25_HEP.RData')
sce
# class: SingleCellExperiment 
# dim: 3073 1362 
# metadata(0):
#   assays(2): counts logcounts
# rownames(3073): Phlda2 Myl7 ... Hmcn1 Tfdp2
# rowData names(2): ENSEMBL SYMBOL
# colnames(1362): cell_63240 cell_63242 ... cell_95715 cell_95717
# colData names(32): cell barcode ... C_Leiden_0.8 C_Leiden_1.2
# reducedDimNames(5): pca.corrected umap TSNE UMAP force
# altExpNames(0):
  

# logmat <- as.matrix(logcounts(sce))
# M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
# save(M, file="CTS_ShrinkM_E8.25_HEP.RData", compress=TRUE) 
# load(file="CTS_ShrinkM_E8.25_HEP.RData")
# dim(M)
# # [1] 3073 3073
# all(rownames(sce) == rownames(M)) #TRUE


######################################################################
##                                                                  ##
## 1) Summarize the CT identifications for                          ##
## the three of the original predicted CT states which are          ## 
## later HEP C15; Endothelium C6, and early HEP S13                 ##
##                                                                  ##
######################################################################

running.methods <- list.dirs()
running.methods <- running.methods[grepl('/C_',running.methods)]
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

## manually add the QuanTC's prediction on softing cluster IDs,
## See S Figures for details
CT$QuanTC.k4_run <- c('C1')
CT$QuanTC.k6_run <- c('C2','C3','C1','C6')


###########################################################################
##                                                                       ##
## 1) Annotate cell identities by comparing to the original cluster IDs. ##
##                                                                       ##
###########################################################################
colnames(colData(sce))
# 
x <- grep("C_", colnames(colData(sce)))
colnames(colData(sce))[x]
 
best.matched.clusterID <- data.frame()
query = c('13','15','6','7') #  the cluster ID for bifurcation which are named 
# as early HEP, later HEP, EP, HP, respectively
# among the studied clusters C3, C13, C10, C6, C15, C7
 
for(i in x){
  tmp <- table(sce$label, colData(sce)[,i])
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
colnames(best.matched.clusterID) <- c('early HEP', 'later HEP', 'EP', 'HP')

## additionally recalculate without the TC among QuanTC's clusters
i= grep('Soft_k4', colnames(colData(sce))) # 7
tmp <- table(sce$label, colData(sce)[,i])
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
                                'C_Soft_k4.wo.TC'= best.match[query])

i= grep('Soft_k6', colnames(colData(sce))) # 7
tmp <- table(sce$label, colData(sce)[,i])
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
                                'C_Soft_k6.wo.TC'= best.match[query])

### add the best matches for the additional BioTIP runs and the original QuanTC run
setdiff(names(CT), rownames(best.matched.clusterID))
# [1] "C_SNNGraph"          "C_SNNGraph_allcells" "QuanTC.k4_run"      "QuanTC.k6_run"

# Now, manually add the best-match infor for these additional BioTIP runs
n <- nrow(best.matched.clusterID) # 14
best.matched.clusterID <- rbind(best.matched.clusterID, 
                                query, query,
                                best.matched.clusterID['C_Soft_k4.wo.TC',],
                                best.matched.clusterID['C_Soft_k6.wo.TC',]) 
rownames(best.matched.clusterID)[(n+1):(n+4)] <- setdiff(names(CT), rownames(best.matched.clusterID))


colnames(best.matched.clusterID) <- c("eHEP", "lHEP", "EP"  ,   "HP")
   
save(best.matched.clusterID, file='best.matched.clusterID.RData')
best.matched.clusterID

################################################################################################
##                                                                                            ##
## 3) calculate Jaccard-similarity score for identifing the 4 potential bufurcations          ##
## among the tested cells of the original cluster IDs C3, C13, C10, C6, C15, C7               ##
## between each predictions and that of our original run of BioTIP with the SNNGraph clusters ##
## after finding the reference cluster ID (i.e. best-matched IDs) for the QuanTC clusters     ##
## And maually adding the 'best matches' for additional BioTIP/QuanTC runs of this datasets   ##
##                                                                                            ##
################################################################################################
n=4  ## 4 potntial bifurcation clusters !!!!!!!!!!!!
jaccard.CT <- array()
for(i in 1:length(CT)){
  x <- names(CT)[i]
  reference <- as.matrix(best.matched.clusterID[x,1:n])
  # # adjust the reference cluster ID for shoft-thresholding clustering
  # if(grepl('C_Soft', x) & !grepl('wo.TC', x) & !grepl('QuanTC_QuanTC', x)) {
  #   reference <- c(reference, best.matched.clusterID[paste0(x,'.wo.TC'),1:n]) %>% unlist()
  # }
  jaccard.CT[i] <- jaccard.sim(unique(reference), CT[[x]]) 
}
names(jaccard.CT) <- names(CT)
round(jaccard.CT, 2)

#############################################################################################
##                                                                                         ##
## 4) we generate a proxy 'golden standard' of CTS genes,                                  ##               
## which were 16 genes identifed by at least 3 out of eight identfications for early HEP   ##                           
##                                                                                         ##
#############################################################################################
{  
  n.select <- c(3,4,4,5)  #c(3,2,2) #  floor(length(select) * 0.8))  
  # for those identified eHEP as a CT cluster
  select <- list()
  select[['eHEP']] <- c( "SNNGraph_allcells", "SNNGraph",    "consensus_ks6",    
                      "Leiden_0.4" , "Leiden_0.8" , "Soft_k4",  "QuanTC.k4_run",    "QuanTC.k6_run")  
  select[['lHEP']] <- c( "SNNGraph_allcells", "consensus_ks4" ,  "consensus_ks8","consensus_ks10",
                         #"consensus_ks5" ,  "consensus_ks7",  "consensus_ks9",        
                         "Leiden_1.2" , "Soft_k4",  "Soft_k4.wo.TC","Soft_k6.wo.TC")  
  select[['EP']] <- c(  "SNNGraph", "consensus_ks8" ,  #  "consensus_ks7", "consensus_ks9",   
                         "Leiden_0.4" , "Leiden_0.8" ,"Leiden_1.2" , "Soft_k6",  "QuanTC.k6_run")  
  select[['HP']] <- c( "SNNGraph_allcells", "SNNGraph", "consensus_ks6" ,    "consensus_ks8",    
                        "Leiden_0.4" , "Leiden_0.8" ,"Soft_k4", "Soft_k4.wo.TC","Soft_k6")  
  
  GS = list()
  for(m in 1:length(select)){
    GS.CTS <- NULL
    for(i in 1:length(select[[m]])){
      if(grepl('QuanTC', select[[m]][i])) {
        load(file=paste0(select[[m]][i],'/CTS.RData'))
        set2 <- CTS  
        GS.CTS = c(GS.CTS, unlist(set2[x]))
        rm(CTS)
      } else  {
        load(file=paste0('C_',select[[m]][i],'/BioTIP.res.RData'))
        set2 <- res$CTS.candidate[which(res$significant)]
        best.match <- best.matched.clusterID[paste0('C_',select[[m]][i]),names(select)[m]]   
        x <- grep(best.match, names(set2))
        GS.CTS = c(GS.CTS, unlist(set2[x]))
      }
    }
    x <- table(GS.CTS)
    GS.CTS <- names(x)[which(x>= n.select[m])]
    GS[[m]] <- GS.CTS
  
    names(GS)[m] <- names(select)[m]
   (n.GS <- length(GS.CTS))  # 40
    write.table(GS[[m]], file=paste0('GS.CTS.',names(select)[m],'_at.least.',n.select[m],'identificaitons_fr.',length(select[[m]]),'predicitons.txt'),
            row.names=F, col.names=F)
  }
  lengths(GS)
 # eHEP       lHEP      EP      HP 
 # 16        64        39       69
  
}


####################################################################
##                                                                ##
## 5) evaluate the F1 score for the identified CTS genes          ##
## for any commonly identifed CTs from any two identifications    ##
##                                                                ##
####################################################################
# load('original_run/CTS.RData')
# lengths(CTS) 
# # 15  6 16 13 # cluster ID
# # 67 90 79 60 # number of CTS genes
# CTS.reference <- CTS[c('13','15',  '6')]
# rm(CTS)
GS.CTS <- read.table('GS.CTS.eHEP_at.least.3identificaitons_fr.8predicitons.txt')
CTS.reference <- GS.CTS[,1]

best.matched.clusterID['C_SNNGraph_allcells',]
#                         early HEP later HEP EP HP
# C_SNNGraph_allcells        13        15   6     7


# load the CTS prediction with different clustering methods
running.methods <- list.dirs()
running.methods <- running.methods[grepl("./C_", running.methods)]
running.methods
# [1] "./C_consensus_ks10"    "./C_consensus_ks4"     "./C_consensus_ks5"    
# [4] "./C_consensus_ks6"     "./C_consensus_ks7"     "./C_consensus_ks8"    
# [7] "./C_consensus_ks9"     "./C_Leiden_0.4"        "./C_Leiden_0.8"       
# [10] "./C_Leiden_1.2"        "./C_SNNGraph"          "./C_SNNGraph_allcells"
# [13] "./C_Soft_k4"           "./C_Soft_k4.wo.TC"     "./C_Soft_k6"          
# [16] "./C_Soft_k6.wo.TC" 
methods <- lapply(running.methods, function(x) unlist(strsplit(x, "./C_"))[2]) %>% unlist()
x <- which(is.na(methods))
if(length(x)>0) { 
  methods <- methods[-x]
  running.methods <- running.methods[-x]
  }
## manually adding the QuanTC's predictions
running.methods <- c(running.methods, './QuanTC.k4_run','./QuanTC.k6_run')
methods <- c(methods,'QuanTC.k4_run','QuanTC.k6_run')

F1.res <- get.multiple.F1.score(match.col.name='eHEP', 
                                CTS.reference = CTS.reference,
                                QuanTC.genes.list.id = 1,
                                running.methods=running.methods, methods=methods,
                                sce = sce,
                                best.matched.clusterID=best.matched.clusterID)
F1.eHEP <- F1.res$F1.scores$eHEP
F1.eHEP.ctl <- F1.res$F1.ctl$eHEP


# # method 1: used, between the inferred GS.CTS and any method
# F1.CTS <- array()
# for(i in 1:length(running.methods)){
#   set1 <- list()
#   set1[['eHEP']]= CTS.reference  
#   x <- methods[i]
#   if(grepl('QuanTC', methods[i] )) {
#     load(file=paste0(running.methods[i],'/CTS.RData'))
#     # only the transition genes for the early HEP were recorded for comparision
#     set2 <- CTS
#     if(length(CTS)>1) set2 <- list('eHEP' = unique(unlist(CTS)))
#     matched.name <- best.matched.clusterID[methods[i], 'eHEP']
#     names(set2) <- names(set1)
#     F1.CTS[i] <- F_score.CTS(set1, set2, weight=TRUE, merge.cluser=TRUE)$F1
#   } else {
#     load(file=paste0(running.methods[i],'/BioTIP.res.RData'))
#     if(length(res)>1){
#       names(res)
#       #[1] "CTS.candidate" "CTS.score"     "Ic.shrink"     "significant" 
#       set2 <- res$CTS.candidate[which(res$significant)]
#       if(length(set2) >0) {
#         matched.name <- best.matched.clusterID[paste0('C_',methods[i]), 'eHEP' ]
#         names(matched.name) <- names(set1)
#         names(set1) <- matched.name[names(set1)]
#         F1.CTS[i] <- F_score.CTS(set1, set2, weight=TRUE)$F1
#       } else  F1.CTS[i] = 0
#     } else  F1.CTS[i] = 0
#   }
# }  
# names(F1.CTS) <- methods
 
round(F1.eHEP,3)
# consensus_ks10     consensus_ks4     consensus_ks5     consensus_ks6     consensus_ks7 
# 0.000             0.000             0.039             0.241             0.000 
# consensus_ks8     consensus_ks9        Leiden_0.4        Leiden_0.8        Leiden_1.2 
# 0.000             0.000             0.341             0.288             0.000 
# SNNGraph SNNGraph_allcells           Soft_k4     Soft_k4.wo.TC           Soft_k6 
# 0.052             0.000             0.044             0.000             0.000 
# Soft_k6.wo.TC     QuanTC.k4_run     QuanTC.k6_run 
# 0.000             0.000             0.000 
# 


## calculating F1 scores for the gene module with the highest DNB score which is C7
# Even for the C13 of the C_SNNGraph run, we can estimate its ability to identify teh robust CTS genes here 
load('C_SNNGraph/BioTIP.res.RData')
lengths(res$CTS.candidate)
#  7  6 13 
# 31 80 45 
F1.eHEP.DNB <- F_score.CTS(list('12'=CTS.reference), res$CTS.candidate[3], weight=TRUE)$F1
F1.eHEP.ctl.DNB <- F_score.CTS(list('2'=CTS.reference), 
                          list('2'=sample(rownames(sce),45)), 
                          weight=TRUE)$F1 

nn = length(F1.eHEP)  # 18
Normalized.F1 <- Normalize.F1(c(F1.eHEP.ctl, F1.eHEP, F1.eHEP.ctl.DNB, F1.eHEP.DNB))
Norm.F1.eHEP.ctl <- Normalized.F1[1:nn] 
Norm.F1.eHEP <- Normalized.F1[(nn+1):(2*nn)] 
Norm.F1.eHEP.ctl.DNB <- Normalized.F1[(2*nn+1):(3*nn)]
Norm.F1.eHEP.DNB <- Normalized.F1[(3*nn+1):(4*nn)] 



###################################################################################################################
##                                                                                                               ##
# 6) reporting table                                                                                             ##
## Verify if the CTS identified the gene of interest which is Etv2                                               ##
## only the 'SNNGraph_allcells' did                                                                                                              ##
###################################################################################################################
# load the CTS prediction with different clustering methods
HEP.marker <-  c('Etv2')   
  
HEP.CTS <- NULL 
for(i in 1:length(running.methods)){
  if(grepl('QuanTC', methods[i] )) {
    load(file=paste0(running.methods[i],'/CTS.RData'))
    set2 <- CTS 
    if(grepl('k4', methods[i])) names(set2) <-'C1' else names(set2) <- rep('C3',2) # normalzie the name in order to use the weight of F1 score
    best.match <- best.matched.clusterID[methods[i],'early HEP']
    x <- grep(best.match, names(set2))
    if(length(x)>0) tmp <- HEP.marker %in% unlist(set2[x]) else tmp <- rep('FALSE', length(HEP.marker))
    HEP.CTS = rbind(HEP.CTS, tmp)
    } else {
      load(file=paste0(running.methods[i],'/BioTIP.res.RData'))
      if(length(res)>1){ # flag 1
        names(res)
        #[1] "CTS.candidate" "CTS.score"     "Ic.shrink"     "significant" 
        set2 <- res$CTS.candidate[which(res$significant)]
      
        best.match <- best.matched.clusterID[paste0('C_',methods[i]),'early HEP']
        x <- grep(best.match, names(set2))
        if(length(x)>0) tmp <- HEP.marker %in% unlist(set2[x]) else tmp <- rep('FALSE', length(HEP.marker))
        HEP.CTS = rbind(HEP.CTS, tmp)
       } else  { # flag 1
        HEP.CTS <- rbind(HEP.CTS, rep('FALSE', length(HEP.marker)))
      }
    }
}
rownames(HEP.CTS) <- methods
colnames(HEP.CTS) <- HEP.marker

save(HEP.CTS, file='HEP.CTS.RData')


##########################################
## reproting table #######################
tb <- data.frame(Jaccard.CT = round(jaccard.CT,3), 
                 F1.eHEP.ctl = round(F1.eHEP.ctl,3), 
                 F1.eHEP = round(F1.eHEP,3),
                 F1.eHEP.ctl.DNB = round(F1.eHEP.ctl.DNB,3), 
                 F1.eHEP.DNB = round(F1.eHEP.DNB,3),
                 Norm.F1.eHEP.ctl = round(Norm.F1.eHEP.ctl,3), 
                 Norm.F1.eHEP = round(Norm.F1.eHEP,3),
                 Norm.F1.eHEP.ctl.DNB = round(Norm.F1.eHEP.ctl.DNB,3), 
                 Norm.F1.eHEP.DNB = round(Norm.F1.eHEP.DNB,3)
                 )
tb <- cbind(tb, HEP.CTS)
nrow(tb) #[1] 18

# reorder to show
myorder <- c('C_SNNGraph_allcells', 'C_SNNGraph',
             paste0('C_consensus_ks',c(4,6,8,10)),
             paste0('C_Leiden_',c(0.4, 0.8, 1.2)),
             'C_Soft_k4','C_Soft_k4.wo.TC', 'C_Soft_k6','C_Soft_k6.wo.TC', 
             'QuanTC.k4_run' ,'QuanTC.k6_run' 
             )
tb$detected.CT <- lapply(CT, toString) %>% unlist()

(mytable <- tb[myorder,
   c('detected.CT','Jaccard.CT', 'F1.eHEP.ctl', 'F1.eHEP',
     'F1.eHEP.ctl.DNB', 'F1.eHEP.DNB',
     'Norm.F1.eHEP.ctl', 'Norm.F1.eHEP',
     'Norm.F1.eHEP.ctl.DNB', 'Norm.F1.eHEP.DNB')])
#                             detected.CT Jaccard.CT F1.CTS F1.nor
# C_SNNGraph_allcells 7, 11, 15, 16, 13, 8      0.429  0.000  0.307
# C_SNNGraph                      7, 6, 13      0.750  0.052  0.485
# C_consensus_ks4                        4      0.333  0.000  0.307
# C_consensus_ks6                     5, 4      0.667  0.241  0.953
# C_consensus_ks8               6, 4, 8, 3      0.400  0.000  0.307
# C_consensus_ks10             6, 2, 10, 4      0.143  0.000  0.307
# C_Leiden_0.4                     6, 2, 5      0.750  0.341  0.995
# C_Leiden_0.8                  7, 6, 9, 4      0.600  0.288  0.982
# C_Leiden_1.2             9, 3, 10, 11, 4      0.286  0.000  0.307
# C_Soft_k4                     TC, C2, C1      0.750  0.044  0.459
# C_Soft_k4.wo.TC                       C2      0.333  0.000  0.307
# C_Soft_k6                     C6, C1, C5      0.400  0.000  0.307
# C_Soft_k6.wo.TC               C1, C5, C4      0.200  0.000  0.307
# QuanTC.k4_run                         C1      0.333  0.000  0.307
# QuanTC.k6_run             C2, C3, C1, C6      0.750  0.000  0.307


write.table(mytable, 
            file='jaccard.CT_F1.CTS.txt', 
            sep='\t')


###############
## PLOT #######
###############

## generate bar plot #######################
df <- read.table('jaccard.CT_F1.CTS.txt', header=T, sep='\t')
tmp <- rownames(df)
tmp = gsub('C_', '', tmp)
df$method = gsub('QuanTBioTIP_run', 'QuanTC_Soft', tmp)
# # show BioTIP_Soft twice for different comarisions (i.e., plot subpannels) and then reorder
# df <- rbind(df, df[c('C_Soft_k4','C_Soft_k6'),] )
# n = nrow(df)
# df[(n-1):n,'method'] <- c("BioTIP.k4_run", "BioTIP.k6_run")
# 
# x <- which(rownames(df) %in% c('C_Soft_k41','C_Soft_k4.wo.TC','C_Soft_k61', 'C_Soft_k6.wo.TC',
#                                'QuanTC.k4_run', 'QuanTC.k6_run', 'C_SNNGraph_allcells'))
# df <- rbind(df[c('C_Soft_k4.wo.TC','C_Soft_k6.wo.TC','C_SNNGraph_allcells'),],
#             df[-x,], df[c('QuanTC.k4_run', 'QuanTC.k6_run', 'C_Soft_k41','C_Soft_k61' ),])
df$method=factor(df$method, levels=df$method)

p1 <- ggbarplot(df, "method", "Jaccard.CT", orientation = "horiz") +
     geom_col(aes(fill = Jaccard.CT)) + 
     scale_fill_gradient2(low = "white",high = "green") +
     theme_bw(base_size=10)  +        # use a black-and-white theme with set font size
    theme(legend.position = "none") + coord_flip()   
p2 <- ggbarplot(df, "method", "Norm.F1.eHEP", orientation = "horiz") +
     geom_col(aes(fill = Norm.F1.eHEP)) + 
     scale_fill_gradient2(low = "white",high = "blue") +
     theme_bw(base_size=10)    +      # use a black-and-white theme with set font size
     theme(legend.position = "none")  + coord_flip()  
# p3 <- ggbarplot(df, "method", "KIT", orientation = "horiz") +
#      geom_col(aes(fill=KIT)) + 
#      scale_fill_manual(values=c("grey", "dark red")) +
#      theme_bw(base_size=10) + 
#      theme(legend.position = "none")   + coord_flip()      

pdf(file='bar_robustness.pdf', width=6, height=4)
gridExtra::grid.arrange(
     p1, p2#, p3  
     ,ncol=2)
dev.off()

