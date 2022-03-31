#### README ############
## There are 4 sections in this code
## Section 1) apply BioTIP to different clusters of 131 cells by different methods             
##          also apply BioTIP to different clusters of 196 cells by cell collection time (age), marker-based clusters, and unsupervised clusters           
## section 2) quantify the robustness         
## section 3) reporting table, compare the results of BioTIP to the resutls of QuanTC  
## section 4) stability of BioTIP's identification
## section 5) TSNE plot for 3.2k HGV, 131 cells along the AT2 trajectory


## last update 2/8/2022
## by Holly Yang  ##############

setwd('F:/projects/BioTIP/result/GSE52583/BioTIP_GSE52583_robustness/')

# the following funciton has been updated on github, but before updating the package to Bioconductor, you need to call it locally to test
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/ForJennifer/optimize.sd_selection.R') 
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/BioTIP.wrap.R')
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/Stability_score.R')
  
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
## focusing on 131 cells along the AT2 trajectory
###################################################
j="AT2.sce.RData"
if(j=="AT2.sce.RData") load('AT2.sce.RData')  else load('sce.RData') #!!!!!!!!!!!!!!!
# format assayNames to run BioTIP
assayNames(sce) <- 'logcounts'

sce
#class: SingleCellExperiment 
# dim: 3198 131 
# metadata(0):
#   assays(1): logcounts
# rownames(3198): 0610007C21Rik 0610007L01Rik ... Zyx l7Rn6
# rowData names(11): ensembl_gene_id gene_short_name ...
# num_cells_expressed use_for_ordering
# colnames(131): E18_1_C08 E18_1_C09 ... E16_1_C94 E16_1_C95
# colData names(17): age cellName ... C_SNNGraph_k8 C_Soft
# reducedDimNames(2): PCA TSNE
# altExpNames(0):

# logmat <- as.matrix(assays(sce)$logcounts)
# M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
# if(j=="AT2.sce.RData") save(M, file="AT2_ShrinkM.RData", compress=TRUE) else save(M, file="ShrinkM.RData", compress=TRUE)
# dim(M) #3198 3198
# rm(logmat)
if(j=="AT2.sce.RData") load(file="AT2_ShrinkM.RData") else load(file="ShrinkM.RData")
dim(M) #3198 3198


####################################################################
## section 1) running BioTIP on different cell clustering outputs ##
####################################################################
{
  ######### 1.1) setting parameters,   ######################
  localHVG.preselect.cut = 0.1 # A positive numeric value setting how many propotion of global HVG to be selected per cell cluster
  getNetwork.cut.fdr = 0.2  # to construct RW network and extract co-expressed gene moduels
  getTopMCI.gene.minsize = 30  # min number of genes in an identified CTS
  getTopMCI.n.states = 2  # A number setting the number of states to check the significance of DNB score (i.e. MCI) 
  # This parameter can be enlarge when inputting more clusters, usually the expected number of transition states +1
  MCIbottom = 2  # A number setting the threshold to prioritize the initially selected CTS candidates.
  # In our experiment, a number between 2 and 4 generated expected resutls.
  
  ############### 1.2) applying BioTIP to all new clustering outputs  ##################
  x <- grep("C_", colnames(colData(sce)))
  colnames(colData(sce))[x]
  # [1] "C_by_marker"     "C_Monocle.k4"    "C_consensus_ks3" "C_consensus_ks5"
  # [5] "C_consensus_ks7" "C_Leiden_0.4"    "C_Leiden_0.8"    "C_Leiden_1.2"   
  # [9] "C_SNNGraph_k10"  "C_SNNGraph_k8"   "C_Soft" 
  
  set.seed(102020)
  for(i in 1:length(x)){
    samplesL <- split(rownames(colData(sce)), f = colData(sce)[x[i]])
    (tmp = lengths(samplesL))
    
     #outputpath = paste0(colnames(colData(sce))[x[i]],'/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
     #load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
     res <- BioTIP.wrap(sce, samplesL, subDir=colnames(colData(sce))[x[i]],
                       getTopMCI.n.states=getTopMCI.n.states, 
                       getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                       MCIbottom=MCIbottom,
                       localHVG.preselect.cut=localHVG.preselect.cut,  #logmat.local.HVG.testres=testres,
                       getNetwork.cut.fdr=getNetwork.cut.fdr, 
                       M=M,
                       permutation.method='both', 
                       empirical.IC.p.cut=0.1, ## when permutating both gene and sampels toestimate the Ic.shrink's significance, 
                       ## we can loose this creterial to get more meanfuling results
                       verbose=TRUE, plot=TRUE)           
    save(res, file=paste0(colnames(colData(sce))[x[i]],'/BioTIP.res.RData'))  
    
  }
  
  ## 1.3) Additional run for soft thresholding without TC
  samplesL <- split(rownames(colData(sce)), f = colData(sce)$C_Soft)
  (tmp = lengths(samplesL))
  # C1 C2 TC 
  # 31 43 57 
  samplesL = samplesL[1:(length(samplesL)-1)]
  
   # outputpath = paste0('C_Soft.wo.TC/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
   # load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
   res <- BioTIP.wrap(sce, samplesL, subDir='C_Soft.wo.TC',
                     getTopMCI.n.states=getTopMCI.n.states, 
                     getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                     MCIbottom=MCIbottom,
                     localHVG.preselect.cut=localHVG.preselect.cut,  #logmat.local.HVG.testres=testres,
                     getNetwork.cut.fdr=getNetwork.cut.fdr, 
                     M=M,
                     permutation.method='both',
                     verbose=TRUE, plot=TRUE)           
  save(res, file='C_Soft.wo.TC/BioTIP.res.RData')  
  
  ## 1.4) additional run for time-course 
  samplesL <- split(rownames(colData(sce)), f = colData(sce)$age)
  (tmp = lengths(samplesL))
  #14.5  16.5  18.5 Adult 
  #45    27    15    44 
  # outputpath = paste0('C_CollectionTime_allcells/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
  # load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
  res <- BioTIP.wrap(sce, samplesL, subDir='C_CollectionTime',
                     getTopMCI.n.states=getTopMCI.n.states, 
                     getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                     MCIbottom=MCIbottom,
                     localHVG.preselect.cut=localHVG.preselect.cut,  #logmat.local.HVG.testres=testres,
                     getNetwork.cut.fdr=getNetwork.cut.fdr, 
                     M=M,
                     permutation.method='both',
                     verbose=TRUE, plot=TRUE)           
  save(res, file='C_CollectionTime/BioTIP.res.RData')  
  
  
  ## 1.5) additional run for all collected time-course cells  
  load('sce.RData')
  sce
  # class: SingleCellExperiment 
  # dim:  3198 196
  samplesL <- split(rownames(colData(sce)), f = colData(sce)$age)
  (tmp = lengths(samplesL))
  # 14.5  16.5  18.5 Adult 
  # 45    27    80    44
  
  # format assayNames to run BioTIP
  assayNames(sce) <- 'logcounts'
  # outputpath = paste0('C_CollectionTime_allcells/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
  # load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
  res <- BioTIP.wrap(sce, samplesL, subDir='C_CollectionTime_allcells',
                     getTopMCI.n.states=getTopMCI.n.states, 
                     getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                     MCIbottom=MCIbottom,
                     ocalHVG.preselect.cut=localHVG.preselect.cut, #logmat.local.HVG.testres=testres,
                     getNetwork.cut.fdr=getNetwork.cut.fdr, 
                     M=M,
                     permutation.method='both',
                     verbose=TRUE, plot=TRUE)           
  save(res, file='C_CollectionTime_allcells/BioTIP.res.RData')  
  
  ## 1.6) additionally run for all collected cells, clustered by gene markers
  samplesL <- split(rownames(colData(sce)), f = colData(sce)$C_by_marker)
  (tmp = lengths(samplesL))
  # Ambiguous       AT1       AT2   Unknown 
  #   6        63        77        50 
  # outputpath = paste0('C_by_marker_allcells/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
  # load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
  res <- BioTIP.wrap(sce, samplesL, subDir='C_by_marker_allcells',
                     getTopMCI.n.states=getTopMCI.n.states, 
                     getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                     MCIbottom=MCIbottom,
                     localHVG.preselect.cut=localHVG.preselect.cut,  # logmat.local.HVG.testres=testres,
                     getNetwork.cut.fdr=getNetwork.cut.fdr, 
                     M=M,
                     permutation.method='both',
                     verbose=TRUE, plot=TRUE)           
  save(res, file='C_by_marker_allcells/BioTIP.res.RData')  
  
  
  ## 1.7) copy QuanTC-identified transition genes
  # BEGIN DO not repeat read and save QuanTC-predicted CTS genes   ##
  CTS <- list()
  tmp <- read.table('../QuanTC_Output/AT2/transition_gene_name1.2.txt',
                    header=FALSE)#E14.5 in soft-thresholing clustering
  CTS[['C1']] <-  tmp[,1]
  
  tmp <- read.table('../QuanTC_Output/AT2/transition_gene_name2.2.txt',
                    header=FALSE)#E16.5 in soft-thresholing clustering
  CTS[['C3']] <- tmp[,1]
  
  save(CTS, file='QuanTC_run/CTS.RData')
  # END DO not repeat read and save QuanTC-predicted CTS genes   ##
  ###################################################################
  
}


####################################################################
## section 2) quantify the robustness                             ##
####################################################################
i= "AT2.sce.RData"
if(j=="AT2.sce.RData") load('AT2.sce.RData')  else load('sce.RData') #!!!!!!!!!!!!!!!
# format assayNames to run BioTIP
assayNames(sce) <- 'logcounts'
{
  ## 2. 1) Summarize the CT identifications
  { 
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
 }
  lapply(CT, toString) %>% unlist()
  # C_by_marker      C_by_marker_allcells          C_CollectionTime 
  # "AT2"                        ""                        "" 
  # C_CollectionTime_allcells           C_consensus_ks3           C_consensus_ks5 
  #         "18.5"                       "1"                       "3" 
  # C_consensus_ks7              C_Leiden_0.4              C_Leiden_0.8 
  #         ""                       "2"                      "NA" 
  # C_Leiden_1.2              C_Monocle.k4            C_SNNGraph_k10 
  #         "NA"                       "2"                        "" 
  # C_SNNGraph_k8                    C_Soft              C_Soft.wo.TC 
  #           "1"                        ""                        "" 
  
  ## manually add the QuanTC's prediction about CT states,
  CT$QuanTC_run <- c('C3','C1') # the soft-thresholding cluster IDs
   
   
  ## 2.2) Annotate cell identities by comparing to expected CT cluster ##
  {
    ## load the QuanTC-identified transition cells
    tmp <- read.table('../QuanTC_Output/AT2/C_TC.txt')
    tmp <- tmp[,1]
    table(tmp)
    #  1  2  3  4  5 
    # 36 43 13  2 37 
    colData(sce)$QuanTC_run = tmp
    ## replace the QuanTC.cluster IDs
    tmp <- data.frame(colData(sce))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 1, 'C1'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 2, 'C2'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 3, 'C3'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 4, 'C4'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 5, 'TC'))
    colData(sce)$QuanTC_run <- tmp$QuanTC_run
    
    colData(sce)$age_cellType <- paste(colData(sce)$age, colData(sce)$putative_cell_type, sep="_")
    table(colData(sce)$QuanTC_run, colData(sce)$age_cellType)
    #    14.5_NA 16.5_NA 18.5_AT2 18.5_BP 18.5_ciliated 18.5_Clara Adult_NA
    # C1      36       0        0       0             0          0        0
    # C2       0       0        0       0             0          0       43
    # C3       0      13        0       0             0          0        0
    # C4       0       0        1       1             0          0        0
    # TC       9      14       10       0             1          2        1
    
    colnames(colData(sce))
  # 
  x <- grep("C_", colnames(colData(sce)))
  colnames(colData(sce))[x]
  
  best.matched.clusterID <- data.frame()
  query = c('16.5_NA','18.5_AT2') #  the age_cellType for bifurcation which are named 
  # as early HEP, later HEP, EP, HP, respectively
  # among the studied clusters C3, C13, C10, C6, C15, C7
  
  for(i in x){
    tmp <- table(sce$age_cellType, colData(sce)[,i])
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
  colnames(best.matched.clusterID) <- c('E16.5', 'E18.5_AT2')
  
  ## additionally recalculate without the TC among QuanTC's clusters
  i= grep('Soft_k4', colnames(colData(sce))) # 7
  tmp <- table(sce$age_cellType, colData(sce)[,i])
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
  # [1] "C_by_marker_allcells"      "C_CollectionTime"         
  # [3] "C_CollectionTime_allcells"
  
  # Now, manually add the best-match infor for these additional BioTIP runs
  n <- nrow(best.matched.clusterID) # 14
  best.matched.clusterID <- rbind(best.matched.clusterID, 
                                  best.matched.clusterID['C_by_marker',],
                                  c('16.5','18.5'),  c('16.5','18.5')) 
  rownames(best.matched.clusterID)[(n+1):(n+3)] <- setdiff(names(CT), rownames(best.matched.clusterID))
  
} 
  save(best.matched.clusterID, file='best.matched.clusterID.RData')
  best.matched.clusterID
  #                             E16.5 E18.5_AT2
  # C_by_marker               Unknown       AT2
  # C_Monocle.k4                    2         1
  # C_consensus_ks3                 1         3
  # C_consensus_ks5                 3         1
  # C_consensus_ks7                 7         2
  # C_Leiden_0.4                    2         1
  # C_Leiden_0.8                    1         5
  # C_Leiden_1.2                    5         6
  # C_SNNGraph_k10                  4         1
  # C_SNNGraph_k8                   1         4
  # C_Soft                         TC        TC
  # QuanTC_run                     C3        TC
  # C_Soft.wo.TC                   C3        C4
  # C_by_marker_allcells      Unknown       AT2
  # C_CollectionTime             16.5      18.5
  # C_CollectionTime_allcells    16.5      18.5
  
 
  ## 2.3) calculate Jaccard-similarity score for identifing the 4 potential bufurcations        ##
  ## among the tested cells of the original cluster IDs C3, C13, C10, C6, C15, C7               ##
  ## between each predictions and that of our original run of BioTIP with the SNNGraph clusters ##
  ## after finding the reference cluster ID (i.e. best-matched IDs) for the QuanTC clusters     ##
  ## And maually adding the 'best matches' for additional BioTIP/QuanTC runs of this datasets   ##
{
  jaccard.CT <- array()
  for(i in 1:length(CT)){
    x <- names(CT)[i]
    n=1  ## 2 potntial bifurcation clusters !!!!!!!!!!!!
    reference <- as.matrix(best.matched.clusterID[x,1:n])
    # # adjust the reference cluster ID for shoft-thresholding clustering
    # if(grepl('C_Soft', x) & !grepl('wo.TC', x) & !grepl('QuanTC_QuanTC', x)) {
    #   reference <- c(reference, best.matched.clusterID[paste0(x,'.wo.TC'),1:n]) %>% unlist()
    # }
    jaccard.CT[i] <- jaccard.sim(unique(reference), CT[[x]]) 
  }
  names(jaccard.CT) <- names(CT)
}
  round(jaccard.CT, 2)
  #  C_by_marker      C_by_marker_allcells          C_CollectionTime 
  # 0                         0                         0 
  # C_CollectionTime_allcells           C_consensus_ks3           C_consensus_ks5 
  # 0                         1                         1 
  # C_consensus_ks7              C_Leiden_0.4              C_Leiden_0.8 
  # 0                         1                         0 
  # C_Leiden_1.2              C_Monocle.k4            C_SNNGraph_k10 
  # 0                         1                         0 
  # C_SNNGraph_k8                    C_Soft              C_Soft.wo.TC 
  # 1                         0                         0 
  # QuanTC_run 
  # 0.5 
  
                                                             
  ## 2.4 we generate a proxy 'golden standard' of CTS genes,                                     
  ## which were identifed by at least 3 out of five identfications for E16.5                            
 ###################################################################################
  {  
 # select <- c('by_marker', 'Leiden_0.4','QuanTC_run','CollectionTime_allcells') #  for E18.5
   n.select <- 3 #  floor(length(select) * 0.8))  
   select <- c( 'consensus_ks3', 'consensus_ks5', 'SNNGraph_k8', 'Leiden_0.4', 'QuanTC_run')  
  
  GS.CTS = NULL
  for(i in 1:length(select)){
    if(select[i] == 'QuanTC_run') {
      load(file=paste0(select[i],'/CTS.RData'))
      set2 <- CTS[[2]]  # CTS[[1]] for E18.5
      GS.CTS = c(GS.CTS, unlist(set2[x]))
      rm(CTS)
    } else  {
      load(file=paste0('C_',select[i],'/BioTIP.res.RData'))
      set2 <- res$CTS.candidate[which(res$significant)]
      best.match <- best.matched.clusterID[paste0('C_',select[i]),'E16.5']  #E18.5_AT2']
      x <- grep(best.match, names(set2))
      GS.CTS = c(GS.CTS, unlist(set2[x]))
    }
  }
  x <- table(GS.CTS)
  GS.CTS <- names(x)[which(x>= n.select)]
  GS.CTS
  # "9130023H24Rik" "Capsl"         "Gm6320"        "Nme5"          "Rsph1"        
  # [6] "Tmem107"  
  (n.GS <- length(GS.CTS))  # 6
  }
  write.table(GS.CTS, file=paste0('GS.CTS.E16.5_at.least.',n.select,'identificaitons_fr.',length(select),'predicitons.txt'),
              row.names=F, col.names=F)
  
  
  ## 2.5) evaluate the F1 score for the identified CTS genes at E16.5    
  ## for any commonly identifed CTs from two identifications    ##
  ################################################################
  CTS.reference <- GS.CTS # inferred reference for this study !!
  {
  # load the CTS prediction with different clustering methods
  running.methods <- list.dirs()
  running.methods <- running.methods[grepl("./C_", running.methods)]
  running.methods
  # [1] "./C_by_marker"               "./C_by_marker_allcells"      "./C_CollectionTime"
  # [4] "./C_CollectionTime_allcells" "./C_consensus_ks3"           "./C_consensus_ks5"
  # [7] "./C_consensus_ks7"           "./C_Leiden_0.4"              "./C_Leiden_0.8"
  # [10] "./C_Leiden_1.2"              "./C_Monocle.k4"              "./C_SNNGraph_k10"
  # [13] "./C_SNNGraph_k8"             "./C_Soft"                    "./C_Soft.wo.TC"

  methods <- lapply(running.methods, function(x) unlist(strsplit(x, "./C_"))[2]) %>% unlist()
  x <- which(is.na(methods))
  if(length(x)>0) {
    methods <- methods[-x]
    running.methods <- running.methods[-x]
  }
  ## manually adding the QuanTC's predictions
  running.methods <- c(running.methods, './QuanTC_run')
  methods <- c(methods,'QuanTC_run')

  F1.res <- get.multiple.F1.score(match.col.name='E16.5',# MCIbottom, 
                                    CTS.reference = GS.CTS,
                                    QuanTC.genes.list.id = 2,
                                    running.methods=running.methods, methods=methods,
                                    sce = sce,
                                    best.matched.clusterID=best.matched.clusterID)
  F1.CTS <- F1.res$F1.scores$E16.5
  F1.CTS.ctl <- F1.res$F1.ctl$E16.5
    
  round(F1.CTS,3)
   # round(F1.CTS,3)
  # by_marker      by_marker_allcells          CollectionTime 
  # 0.000                   0.000                   0.000 
  # CollectionTime_allcells           consensus_ks3           consensus_ks5 
  # 0.000                   0.000                   0.162 
  # consensus_ks7              Leiden_0.4              Leiden_0.8 
  # 0.000                   0.136                   0.000 
  # Leiden_1.2              Monocle.k4            SNNGraph_k10 
  # 0.000                   0.000                   0.000 
  # SNNGraph_k8                    Soft              Soft.wo.TC 
  # 0.176                   0.000                   0.000 
  # QuanTC_run 
  # 0.000 
  

  ## calculating F1 scores for the gene module with the highest DNB score
  # which is C2_38g of the Leiden_0.4 run
  load('C_Leiden_0.4/BioTIP.res.RData')
  lengths(res$CTS.candidate)
  # 2  1  2  1 
  # 38 67 44 39 
  intersect(GS.CTS, res$CTS.candidate[[1]])
  #character(0)
  F1.DNB <- F_score.CTS(list('2'=GS.CTS), res$CTS.candidate[1], weight=TRUE)$F1
  F1.ctl.DNB <- F_score.CTS(list('2'=GS.CTS), 
                            list('2'=sample(rownames(sce),38)), 
                            weight=TRUE)$F1 
  
  nn = length(F1.CTS)  # 16
  Normalized.F1 <- Normalize.F1(c(F1.CTS.ctl, F1.CTS, F1.ctl.DNB, F1.DNB))
  Norm.F1.CTS.ctl <- Normalized.F1[1:nn] 
  Norm.F1.CTS <- Normalized.F1[(nn+1):(2*nn)] 
  Norm.F1.ctl.DNB <- Normalized.F1[(2*nn+1)]
  Norm.F1.DNB <- Normalized.F1[(2*nn+2)] 
  
}

###########################################################
## Section 3) reproting table, compare BioTIP to QuanTC  ##
###########################################################
{
    tb <- data.frame(Jaccard.CT = round(jaccard.CT,3), 
                     F1.CTS.ctl = round(F1.CTS.ctl,3), 
                     F1.CTS = round(F1.CTS,3), 
                     Norm.F1.CTS.ctl = round(Norm.F1.CTS.ctl,3), 
                     Norm.F1.CTS = round(Norm.F1.CTS,3), 
                     F1.ctl.DNB = round(F1.ctl.DNB,3),
                     F1.DNB = round(F1.DNB,3),
                     Norm.F1.ctl.DNB = round(Norm.F1.ctl.DNB,3),
                     Norm.F1.DNB = round(Norm.F1.DNB,3)
                     )
  nrow(tb) #[1] 16
  
  # reorder to show
  myorder <- c('C_by_marker', 'C_CollectionTime',
               paste0('C_consensus_ks',c(3,5,7)),
               paste0('C_Leiden_',c(0.4, 0.8, 1.2)),
               'C_SNNGraph_k8','C_SNNGraph_k10', 
               'QuanTC_run' 
  )
  tb$detected.CT <- lapply(CT, toString) %>% unlist()
  
  (mytable <- tb[myorder,
                 c('detected.CT','Jaccard.CT', 
                   'F1.CTS.ctl', 'F1.CTS', 'F1.ctl.DNB','F1.DNB',
                   'Norm.F1.CTS.ctl', 'Norm.F1.CTS', 'Norm.F1.ctl.DNB','Norm.F1.DNB')])
  #                 detected.CT Jaccard.CT F1.CTS F1.nor
  # C_by_marker              AT2          0  0.000  0.322
  # C_CollectionTime                      0  0.000  0.322
  # C_consensus_ks3            1          1  0.000  0.322
  # C_consensus_ks5            3          1  0.162  0.980
  # C_consensus_ks7                       0  0.000  0.322
  # C_Leiden_0.4               2          1  0.136  0.952
  # C_Leiden_0.8              NA          0  0.000  0.322
  # C_Leiden_1.2              NA          0  0.000  0.322
  # C_SNNGraph_k8              1          1  0.176  0.989
  # C_SNNGraph_k10                        0  0.000  0.322
  # QuanTC_run                C3          1  0.000  0.322
  
  
 ## in fact, we only estimated for DNB using Leiden_0.4 clustering 
  mytable[which(rownames(mytable) !='C_Leiden_0.4'),
          c('Norm.F1.ctl.DNB', 'Norm.F1.DNB')]  <- NA
  write.table(mytable, 
              file='jaccard.CT_F1.CTS.txt', 
              sep='\t')
  
  
  ###############
  ## PLOT #######
  ###############
  
  ## generate bar plot #######################
  df <- read.table('jaccard.CT_F1.CTS.txt', header=T, sep='\t')
  df <- df[, -which(colnames(df) %in% c('Norm.F1.ctl.DNB', 'Norm.F1.DNB'))]
  tmp <- rownames(df)
  tmp = gsub('C_', '', tmp)
  df$method=tmp
  
  df$method=factor(df$method, levels=df$method)
  
  p1 <- ggbarplot(df, "method", "Jaccard.CT", orientation = "horiz") +
    geom_col(aes(fill = Jaccard.CT)) + 
    scale_fill_gradient2(low = "white",high = "green") +
    theme_bw(base_size=10)  +        # use a black-and-white theme with set font size
    theme(legend.position = "none") + coord_flip() 
  
  p2 <- ggbarplot(df, "method", "Norm.F1.CTS", orientation = "horiz") +
    geom_col(aes(fill = Norm.F1.CTS)) + 
    scale_fill_gradient2(low = "white",high = "blue") +
    theme_bw(base_size=10)    +      # use a black-and-white theme with set font size
    theme(legend.position = "none")  + coord_flip()  
   
  pdf(file='bar_robustness.pdf', width=6, height=4)
  gridExtra::grid.arrange(
    p1, p2#, p3  
    ,ncol=2)
  dev.off()
  }
  
  
###########################################################
## Section 4) stability using the Leiden_0.4 clusters    ##
###########################################################
{
  if(j=="AT2.sce.RData") load(file="AT2_ShrinkM.RData") else load(file="ShrinkM.RData")
dim(M) #3198 3198

colnames(colData(sce))
# [1] "age"                "cellName"           "GEO.access"         "putative_cell_type"
# [5] "C_by_marker"        "C_Monocle.k4"       "Pseudotime"         "Size_Factor"       
# [9] "C_consensus_ks3"    "C_consensus_ks5"    "C_consensus_ks7"    "C_Leiden_0.4"      
# [13] "C_Leiden_0.8"       "C_Leiden_1.2"       "C_SNNGraph_k10"     "C_SNNGraph_k8"     
# [17] "C_Soft"   
sce$label <- sce$C_Leiden_0.4  # selected clustering method
rowData(sce) <- NULL # to allow run sample()

######### 4.1) setting parameters, the same as section 1  ######################
## but try multiple getTopMCI.gene.minsize 
getTopMCI.gene.minsize <- c(30,20,10)

######### 4.2) run BioTIP on downsized data ##################
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
      this.run <- BioTIP.wrap(sce.dnsize, samplesL, subDir='stability', globa.HVG.select=FALSE,
                              getTopMCI.n.states=getTopMCI.n.states, getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                              n.getMaxMCImember = 2, MCIbottom=MCIbottom,
                              local.HVG.optimize =TRUE, localHVG.preselect.cut=localHVG.preselect.cut, localHVG.runs=100, #logmat.local.HVG.testres=testres,
                              getNetwork.cut.fdr=getNetwork.cut.fdr, 
                              n.permutation = 100,  M=M,
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
    save(res, flag, counts, file=paste0('stability/BioTIP_res_20runs_MCIbottom',MCIbottom[i],
                                        '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    
  }
}

## 4.2) evaluate the stability of CT detection
## using the robustly detected C2 as reference
###########################################################
df.CT =NULL
for(m in 1:length(getTopMCI.gene.minsize)){
  CT.list <- list()
  Freq <- matrix(0, nrow=nlevels(sce$label), ncol=length(MCIbottom))#!!!
  rownames(Freq) <- levels(sce$label) #!!!
  colnames(Freq) <- MCIbottom 
  for(i in 1:length(MCIbottom)){
    load(file=paste0('stability/BioTIP_res_20runs_MCIbottom',MCIbottom[i],
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
  # save(CT.list, file= paste0('stability/CT.detection_variable.MCIbottom_Modulesize',getTopMCI.gene.minsize[m],'.Rdata'))
}
df.CT$Frequency = df.CT$Frequency/n.scability
save(CT.list, df.CT, 
     file= 'stability/CT.detection_variable.MCIbottom_Modulesize.Rdata')


## 4.3) plot for variable Modulesize, when MCIbottom==2         ##
#############################################################
load(file= 'stability/CT.detection_variable.MCIbottom_Modulesize.Rdata')
df.CT$moduleSize <- as.character(df.CT$moduleSize )
ggplot(data=subset(df.CT, MCIbottom==2),
       aes(x=clusterID, y=Frequency, group=moduleSize)) +
  geom_line(aes(color=moduleSize), size=1.5)+
  geom_point(aes(color=moduleSize)) 
dev.copy2pdf(file='stability/CT.detection_variable.MCIbottom2_Modulesize.pdf')


## 4.4) evaluate the stability of CTS identification
## for the C2
###########################################################
load(file= 'stability/CT.detection_variable.MCIbottom_Modulesize.Rdata')
names(CT.list)
#[1] "2"  

# load('original_run/CTS.RData')
# names(CTS)
# #[1] "9"  "10" "9" 
# C9.ref = CTS[which(names(CTS)=='9')]
# C10.ref = CTS['10']
GS.CTS <- read.table('GS.CTS.E16.5_at.least.3identificaitons_fr.5predicitons.txt')
C2.ref <- list('2'=GS.CTS[,1])

MCIbottom  # 2
i=1
F1.C2 <- F1.C2.ctl <- list()
for(m in 1:length(getTopMCI.gene.minsize)){
  load(file=paste0('stability/BioTIP_res_20runs_MCIbottom',MCIbottom[i],
                   '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
  F1.CTS.C2 <- F1.bk.C2 <- array(dim=n.scability)
  for(j in 1:n.scability){
    # evaluate the F1 score for two CTS identificitons
    set2 <- res[[j]]$CTS.candidate[which(res[[j]]$significant)] 
    F1.CTS.C2[j] <- F_score.CTS(C2.ref, set2, weight=TRUE)$F1
    # negative control
    n.module <- lengths(set2)
    random.set2 <- lapply(n.module, function(x) sample(rownames(sce), x))
    F1.bk.C2[j] <- F_score.CTS(C2.ref, random.set2, weight=TRUE)$F1
  } 
  F1.C2[[m]] =  F1.CTS.C2
  F1.C2.ctl[[m]] =  F1.bk.C2
 }

names(F1.C2) <- names(F1.C2.ctl) <- getTopMCI.gene.minsize

save(F1.C2, F1.C2.ctl,
     file='stability/F1.CTS_MCIbottom2.RData')

lapply(F1.C2, mean) %>% unlist()
#        10         20         30 
# 0.12294563 0.12294563 0.08011003
lapply(F1.C2.ctl, mean) %>% unlist()
#0.0024092971 0.0039024701 0.0008333333 

df <- data.frame(F1 = c(unlist(F1.C2.ctl),unlist(F1.C2)),
                 minModuleSize = c(rep(getTopMCI.gene.minsize, lengths(F1.C2.ctl)),
                               rep(getTopMCI.gene.minsize, lengths(F1.C2))),
                 CT = c(rep('C2', sum(lengths(F1.C2.ctl))+sum(lengths(F1.C2)))),
                 CTS = c(rep('ctl', sum(lengths(F1.C2.ctl))), rep('CTS', sum(lengths(F1.C2)))) 
)
df$minModuleSize <- as.character(df$minModuleSize) 
df$CTS <- factor(df$CTS, levels = c('ctl','CTS')) 
df$color.lab <- paste0(df$minModuleSize,df$CTS)
df$color.lab <- factor(df$color.lab, levels=c('10ctl','10CTS','20ctl','20CTS','30ctl','30CTS'))


# To ensure that easy and difficult CT state have equal influence on the final score, 
# we normalized the scores at each CT state across the different MCIbottom setting.
df$Norm.F1 <- Normalize.F1(df$F1)

library(ggbeeswarm)
library(grid)
library(gridExtra)
library(ggpubr)


p1 <- ggplot(data=df, aes(x=CT, y=F1, color=color.lab)) +
  geom_violin(position = position_dodge(width = 0.5)) +
  geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
  xlab('predicted for cluster') + ylab('F1 score for the CTS') +
  scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue")) 

p2 <- ggplot(data=df, aes(x=CT, y=Norm.F1, color=color.lab)) +
  geom_violin(position = position_dodge(width = 0.5)) +
  geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
  xlab('predicted for cluster') + ylab('Normalized F1 score for the CTS') +
  scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue")) 

pdf(file='stability/CTS.F1_variable.minModuleSize.pdf')
gridExtra::grid.arrange(
  p1, p2, top = "CTS stability", bottom = "different module size settings"
  ,ncol=1)
dev.off()

tmp <- compare_means(Norm.F1 ~ color.lab, data=df, group='CT', method='t.test')
subset(tmp, group1=='10ctl' & group2=='10CTS')
# # A tibble: 2 x 9
## A tibble: 1 x 9
# CT    .y.     group1 group2          p    p.adj p.format p.signif method
# <chr> <chr>   <chr>  <chr>       <dbl>    <dbl> <chr>    <chr>    <chr> 
#   1 C2    Norm.F1 10ctl  10CTS  0.00000240 0.000031 2.4e-06  ****     T-test
subset(tmp, group1=='20ctl' & group2=='20CTS')
# A tibble: 1 x 9
# CT    .y.     group1 group2          p    p.adj p.format p.signif method
# <chr> <chr>   <chr>  <chr>       <dbl>    <dbl> <chr>    <chr>    <chr> 
#   1 C2    Norm.F1 20ctl  20CTS  0.00000285 0.000031 2.8e-06  ****     T-test
subset(tmp, group1=='30ctl' & group2=='30CTS')
# A tibble: 1 x 9
# CT    .y.     group1 group2        p  p.adj p.format p.signif method
# <chr> <chr>   <chr>  <chr>     <dbl>  <dbl> <chr>    <chr>    <chr> 
#   1 C2    Norm.F1 30ctl  30CTS  0.000284 0.0026 0.00028  ***      T-test

}
###################################################################
## section 5) TSNE plot for 3.2k HGV, 131 cells along the AT2 trajectory
###################################################################
{
  fid <- "C_Leiden_0.4/"
  load(file=paste0(fid,"BioTIP_top0.1FDR0.2_optimized_local.HVG_selection.RData"))
  names(testres)
  #[1] "1" "2" "3"
  
  library(scater)
  pdf(file='tsne_131cells_C_Leiden_0.4.pdf')
  plotReducedDim(sce, dimred='TSNE', shape_by='age', 
                 colour_by = "C_Leiden_0.4", 
                 text_by='putative_cell_type', 
                 text_size = 4, text_colour='black', point_size=2) + 
    ggtitle('putative_cell_type')  #+ ylim(5,20)  
  
  plotReducedDim(sce, dimred='TSNE', shape_by='age', 
                 colour_by = 'putative_cell_type', 
                 text_by="C_Leiden_0.4", 
                 text_size = 4, text_colour='black', point_size=2) + 
    ggtitle('putative_cell_type')  #+ ylim(5,20)  
  
  
  dev.off()
  
}


