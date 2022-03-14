#### README ############
## There are 4 sections in this code
## Section 1) apply BioTIP to predefiend cellsubtype and cell type, respectively             
## section 2) quantify the robustness after running BioTIP on cells clustered by different methods (parameters)     
## section 3) compare the results of BioTIP to the resutls of QuanTC  
## section 4) stability of BioTIP's identification by iteratively running BioTIP on the downsampled datasets   


## last update 2/17/2022
## by Holly Yang  ##############

#setwd("F:/projects/scRNA/results/AJ/GSE130146_xy/Results_1k/LibrarySize/GSE130146_robustness/") 
setwd("E:/Git_Holly/scRNAseq_examples/result/EB_Zhao2019")

# the following funciton has been updated on github, but before updating the package to Bioconductor, you need to call it locally to test
#source(file='F:/projects/scRNA/source/GSE87038_gastrulation/ForJennifer/optimize.sd_selection.R') 
source(file='E:/Git_Holly/scRNAseq_examples/code/BioTIP.wrap.R')
source(file='E:/Git_Holly/scRNAseq_examples/code/Stability_score.R')
  
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
###################################################
load('../../data/EB_Zhao2019/sce.GSE130146_noenderdormPgerm.RData')
sce
# class: SingleCellExperiment 
# dim: 4000 1531 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(4000): ENSMUSG00000030544 ENSMUSG00000060461 ...
# ENSMUSG00000070814 ENSMUSG00000031951
# rowData names(2): ID Symbol
# colnames(1531): AAACCTGAGTTTCCTT-1 AAACCTGCATGAAGTA-1 ...
# TTTGTCATCGGCTTGG-1 TTTGTCATCTGCTTGC-1
# colData names(13): Sample Barcode ... C_Leiden_0.8 C_Leiden_1.2
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):
 
rownames(sce) <- rowData(sce)$Symbol
rownames(sce) [1:4]
#[1] "Mesp1"  "Dppa5a" "Phlda2" "Tdgf1" 

# logmat <- as.matrix(assays(sce)$logcounts)
# M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
# save(M, file="ShrinkM.RData", compress=TRUE)
# dim(M) #4000 4000
# rm(logmat)
load(file="ShrinkM.RData")
dim(M) #4000 4000



####################################################################
## section 1) running BioTIP on different cell clustering outputs ##
####################################################################
{
  ######### 1.1) setting parameters,   ######################
  localHVG.preselect.cut = 0.1 # A positive numeric value setting how many propotion of global HVG to be selected per cell cluster
  getNetwork.cut.fdr = 0.2  # to construct RW network and extract co-expressed gene moduels
  getTopMCI.gene.minsize = 30  # min number of genes in an identified CTS
  getTopMCI.n.states = 9  # A number setting the number of states to check the significance of DNB score (i.e. MCI) 
  # This parameter can be enlarge when inputting more clusters, usually the expected number of transition states +1
  MCIbottom = c(1.5, 2)  # A number setting the threshold to prioritize the initially selected CTS candidates.
  # In our experiment, a number between 2 and 4 generated expected resutls.
  
  ############### 1.2) applying BioTIP to all new clustering outputs  ##################
  x <- grep("C_", colnames(colData(sce)))
  colnames(colData(sce))[x]
  # [1] "C_SNNGraph_k5"    "C_SNNGraph_k10"   "C_consensus_ks13" "C_consensus_ks15"
  # [5] "C_consensus_ks17" "C_consensus_ks19"  "C_Leiden_0.4"     "C_Leiden_0.8"    
  # [9] "C_Leiden_1.2"   "C_Soft"
  
  
  set.seed(102020)
  for(i in 1:length(x)){
    samplesL <- split(rownames(colData(sce)), f = colData(sce)[x[i]])
    (tmp = lengths(samplesL))
    
    # outputpath = paste0(colnames(colData(sce))[x[i]],'/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
    # load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
    res <- BioTIP.wrap(sce, samplesL, subDir=colnames(colData(sce))[x[i]],
                       getTopMCI.n.states=getTopMCI.n.states, 
                       getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                       MCIbottom=MCIbottom,
                       localHVG.preselect.cut=localHVG.preselect.cut, 
                       # logmat.local.HVG.testres=testres,
                       getNetwork.cut.fdr=getNetwork.cut.fdr, 
                       M=M,
                       verbose=TRUE, plot=TRUE)           
    save(res, file=paste0(colnames(colData(sce))[x[i]],'/BioTIP.res.RData'))  
    
  }
  
  ############### for k=8, soft clustering generate 9 clusters  ##################
  samplesL <- split(rownames(colData(sce)), f = colData(sce)$C_Soft)
  (tmp = lengths(samplesL))
  #  C1  C2  C3  C4  C5  C6  C7  C8  TC 
  # 230 154 157 159 126  92  72  99 442 
  ## 2nd run: without TC 
  samplesL = samplesL[1:8]
  # outputpath = paste0(colnames(colData(sce))[x[i]],'/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
  # load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
  res <- BioTIP.wrap(sce, samplesL, subDir='C_Soft_wo.TC',
                     getTopMCI.n.states=getTopMCI.n.states, 
                     getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                     MCIbottom=MCIbottom,
                     localHVG.preselect.cut=localHVG.preselect.cut, 
                     # logmat.local.HVG.testres=testres,
                     getNetwork.cut.fdr=getNetwork.cut.fdr, 
                     M=M,
                     verbose=FALSE, plot=TRUE)           
  
  save(res, file='C_Soft_wo.TC/BioTIP.res.RData')  
  
 
  ## 1.3) copy QuanTC-identified transition genes, the first two CT# were round the FLK1+ population
  # BEGIN DO not repeat read and save QuanTC-predicted CTS genes   ##
  CTS <- list()
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name1.1.txt',header=FALSE)
  CTS[['C1-C4']] <-  tmp[,1]
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name1.1.txt', header=FALSE)
  CTS[['C1-C4']] <-  c(CTS[['C1-C4']], tmp[,1])
  
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name2.1.txt', header=FALSE)
  CTS[['C4-C2']] <- tmp[,1]
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name2.2.txt',header=FALSE)
  CTS[['C4-C2']] <- c(CTS[['C4-C2']], tmp[,1])
  
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name4.1.txt', header=FALSE)
  CTS[['C6-C5']] <- tmp[,1]
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name4.2.txt',header=FALSE)
  CTS[['C6-C5']] <- c(CTS[['C6-C5']], tmp[,1])
  
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name5.1.txt', header=FALSE)
  CTS[['C5-C8']] <- tmp[,1]
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name5.2.txt',header=FALSE)
  CTS[['C5-C8']] <- c(CTS[['C5-C8']], tmp[,1])
  
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name6.1.txt', header=FALSE)
  CTS[['C8-C3']] <- tmp[,1]
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name6.2.txt',header=FALSE)
  CTS[['C8-C3']] <- c(CTS[['C8-C3']], tmp[,1])
  
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name7.1.txt', header=FALSE)
  CTS[['C7']] <- tmp[,1]
  tmp <- read.table('../QuanTC_Output/k8/transition_gene_name7.2.txt',header=FALSE)
  CTS[['C7']] <- c(CTS[['C7']], tmp[,1])
  
  save(CTS, file='QuanTC_run/CTS.RData')
  # END DO not repeat read and save QuanTC-predicted CTS genes   ##
  ###################################################################
  
}


####################################################################
## section 2) quantify the robustness                             ##
####################################################################
load('../../data/EB_Zhao2019/sce.GSE130146_noenderdormPgerm.RData')
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
  #             C_consensus_ks13                 C_consensus_ks15 
  # "10, 5, 12, 6, 4, 7, 13, 9" "7, 14, 12, 4, 15, 11, 6, 13, 8" 
  # C_consensus_ks17                 C_consensus_ks19 
  # "10, 3, 7, 14, 8, 16, 2"        "6, 5, 16, 10, 19, 8, 14" 
  # C_Leiden_0.4                     C_Leiden_0.8 
  # "6, 5, 4"       "12, 4, 10, 9, 8, 3, 6, 1" 
  # C_Leiden_1.2                   C_SNNGraph_k10 
  # "13, 10, 3, 11, 8, 6, 9, 7"               "9, 6, 2, 1, 4, 3" 
  # C_SNNGraph_k5                           C_Soft 
  # "5, 3, 6, 10, 11, 4, 7"                         "C8, C5" 
  # C_Soft_wo.TC 
  # "C7, C8, C5, C6" 
  
  
  ## manually add the QuanTC's prediction about CT states,
   CT$QuanTC_run <- c('C1','C4', 'C2', 'C6', 'C5', 'C8', 'C3', 'C7')
   
   
  # ## 2.2) Annotate cell identities by comparing to expected CT cluster ##
  {
    ## load the QuanTC-identified transition cells
    tmp <- read.table('../QuanTC_Output/k8/C_TC.txt')
    tmp <- tmp[,1]
    table(tmp)
    #   1   2   3   4   5   6   7   8   9 
    # 230 154 157 159 126  92  72  99 442
    colData(sce)$QuanTC_run = tmp
    ## replace the QuanTC.cluster IDs
    tmp <- data.frame(colData(sce))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 1, 'C1'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 2, 'C2'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 3, 'C3'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 4, 'C4'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 5, 'C5'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 6, 'C6'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 7, 'C7'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 8, 'C8'))
    tmp <- tmp %>% mutate(QuanTC_run = replace(QuanTC_run, QuanTC_run == 9, 'TC'))
    colData(sce)$QuanTC_run <- tmp$QuanTC_run

    colnames(colData(sce))
  
  x <- grep("C_", colnames(colData(sce)))
  colnames(colData(sce))[x]
  # [1] "C_SNNGraph_k5"    "C_SNNGraph_k10"   "C_consensus_ks13" "C_consensus_ks15"
  # [5] "C_consensus_ks17" "C_consensus_ks19" "C_Leiden_0.4"     "C_Leiden_0.8"    
  # [9] "C_Leiden_1.2"     "C_Soft"           "QuanTC_run" 
 
  best.matched.clusterID <- data.frame()
  query = c('11','6','4') #  refere to C_SNNGraph_k5
  # as mesoderm-endoderm progenitor,	early multipotential,	cardiac progenitor, respectively
  # among the studied 13 clusters
  
  for(i in x){
    tmp <- table(sce$C_SNNGraph_k5, colData(sce)[,i])
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
  colnames(best.matched.clusterID) <- c('MEP',	'MP',	'CMP') 
  
  ## additionally recalculate without the TC among QuanTC's clusters
  i= grep('Soft', colnames(colData(sce))) # 7
  tmp <- table(sce$C_SNNGraph_k5, colData(sce)[,i])
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
                                  'C_Soft_wo.TC'= best.match[query])
 
  ### add the best matches for the additional BioTIP runs and the original QuanTC run
  setdiff(names(CT), rownames(best.matched.clusterID))
  #NA
  
  n <- nrow(best.matched.clusterID) # 14

} 
  save(best.matched.clusterID, file='best.matched.clusterID.RData')
  best.matched.clusterID
  #                     MEP MP CMP
  # C_SNNGraph_k5       11    6    4
  # C_SNNGraph_k10       3    2    8
  # C_consensus_ks13     4    6    9
  # C_consensus_ks15    11    7    8
  # C_consensus_ks17    11    7   15
  # C_consensus_ks19    11    6   15
  # C_Leiden_0.4         2    4    3
  # C_Leiden_0.8         1    3    6
  # C_Leiden_1.2         1    7   13
  # C_Soft              TC   C3   C5
  # QuanTC_run          TC   C3   C5
  # C_Soft_wo.TC        C2   C3   C5
 
  ## 2.3) calculate Jaccard-similarity score for identifing the 3 potential bufurcations        ##
  ## among the tested cells of the cluster IDs of "C_SNNGraph_k5" method             ##
  ## between each predictions and that of the run of BioTIP with the "C_SNNGraph_k5" clusters ##
  ## after finding the reference cluster ID (i.e. best-matched IDs) for the QuanTC clusters     ##
  ## And maually adding the 'best matches' for additional BioTIP/QuanTC runs of this datasets   ##
{
  jaccard.CT <- array()
  for(i in 1:length(CT)){
    x <- names(CT)[i]
    n=3  ## 2 potntial bifurcation clusters !!!!!!!!!!!!
    reference <- as.matrix(best.matched.clusterID[x,1:n])
    jaccard.CT[i] <- jaccard.sim(unique(reference), CT[[x]]) 
  }
  names(jaccard.CT) <- names(CT)
}
  round(jaccard.CT, 2)
  # C_consensus_ks13 C_consensus_ks15 C_consensus_ks17 C_consensus_ks19 
  # 0.38             0.33             0.11             0.11 
  # C_Leiden_0.4     C_Leiden_0.8     C_Leiden_1.2   C_SNNGraph_k10 
  # 0.20             0.38             0.22             0.29 
  # C_SNNGraph_k5           C_Soft     C_Soft_wo.TC       QuanTC_run 
  # 0.43             0.25             0.17             0.22 
  
  
                                                             
  ## 2.4.1 we generate a proxy 'golden standard' of CTS genes,                                     
  ## which were identifed by at least 3 out of five identfications
  ## for mesoderm-endoderm bifurcation  (C11 in C_SNNGraph_k5) -----------------------------------
 
    
  load(file='best.matched.clusterID.RData')
  n.select <- 3 #  floor(length(select) * 0.8))  
  select <- c( 'consensus_ks13', 'consensus_ks15', 'SNNGraph_k5','SNNGraph_k10','Leiden_0.8')  
  GS.CTS <- find.proxy.GS(n.select, select,best.matched.clusterID, 'MEP')
  GS.CTS
  # [1] "1700019B21Rik" "Arhgap36"      "Bmp2"          "Cd40"         
  # [5] "Cdkn2b"        "Cnrip1"        "Coch"          "Crym"         
  # [9] "Dhrs7"         "Dusp2"         "Etv2"          "Exoc3l"       
  # [13] "Eya1"          "Fam84b"        "Frmd4b"        "Gata2"        
  # [17] "Gm2694"        "Gng11"         "Lmo2"          "Mmp9"         
  # [21] "Myc"           "Myzap"         "Nav1"          "Nkx2-5"       
  # [25] "Nxph4"         "Phox2a"        "Pwwp2b"        "Rhoj"         
  # [29] "Rnase4"        "Sapcd1"        "Slc39a4"       "Stat3"        
  # [33] "Tal1"          "Tbc1d14"       "Trim11"        "Vax1" 
  (n.GS <- length(GS.CTS))  # 36
  write.table(GS.CTS, file=paste0('GS.CTS.MEP_at.least.',n.select,'identificaitons_fr.',length(select),'predicitons.txt'),
              row.names=F, col.names=F)
  

  ## 2.4.2 we generate a proxy 'golden standard' of CTS genes,                                     
  ## which were identifed by at least 4 out of nine identfications
  ## for the early pluriprotential bifurcation  --------------------------
    n.select <- 4 #  floor(length(select) * 0.8))  
    select <- c("SNNGraph_k5" ,   "SNNGraph_k10",   "consensus_ks13", "consensus_ks15",
                 "consensus_ks17", "Leiden_0.4" ,    "Leiden_0.8" ,  "Leiden_1.2",
               "QuanTC_run"  # the same resutls when including or excluding this 
                )  
    GS.CTS <- find.proxy.GS(n.select, select,best.matched.clusterID, "MP")
    GS.CTS
    # [1] "1190005I06Rik" "2410141K09Rik" "Abcc4"         "Apln"         
    # [5] "Arhgap8"       "Arhgef16"      "Atp1b1"        "BC033916"     
    # [9] "Bspry"         "Calca"         "Cbr3"          "Cbx7"         
    # [13] "Crb3"          "Dppa2"         "Dppa4"         "Enpp3"        
    # [17] "Fgfbp1"        "Foxd3"         "Gdf3"          "Gm13151"      
    # [21] "Gm13154"       "Gm16233"       "Gng3"          "Gtsf1l"       
    # [25] "Ltbp4"         "Ooep"          "Pcdh1"         "Pdk1"         
    # [29] "Pou3f1"        "Ppp1r1a"       "Rbmxl2"        "RP23-362B7.1" 
    # [33] "RP24-105B23.1" "Slc7a3"        "Sox3"          "Syngr1"       
    # [37] "Tmem54"        "Trim61"        "Tubb3"         "Uchl1"        
    # [41] "Utf1"          "Zfp462"        "Zfp534"  
    (n.GS <- length(GS.CTS))  # 43
  write.table(GS.CTS, file=paste0('GS.CTS.MP_at.least.',n.select,'identificaitons_fr.',length(select),'predicitons.txt'),
              row.names=F, col.names=F)
  
  ## 2.4.2 we generate a proxy 'golden standard' of CTS genes,                                     
  ## which were identifed by at least 3 out of 6 identfications
  ## for the early pluriprotential bifurcation  --------------------------
  n.select <- 3 #  floor(length(select) * 0.8))  
  select <- c("SNNGraph_k5" ,    "consensus_ks13", "consensus_ks15",
               "Leiden_0.8" ,  "Leiden_1.2",
              "QuanTC_run"  # the same resutls when including or excluding this 
  )  
  GS.CTS <- find.proxy.GS(n.select, select,best.matched.clusterID, "CMP")
  GS.CTS
  # [1] "Amph"    "Apobec2" "Bnipl"   "Cd74"    "Clu"     "Cpm"     "Dnaaf3" 
  # [8] "Elf3"    "Emp2"    "Gbp3"    "Gbp4"    "Ifi27"   "Ifi35"   "Igsf3"  
  # [15] "Il17rc"  "Irgm1"   "Isg15"   "Kdelr3"  "Mfap2"   "Mov10"   "Nmi"    
  # [22] "Nptx2"   "Parp12"  "Parp9"   "Pdgfrl"  "Pgam2"   "Psmb9"   "Rab20"  
  # [29] "Rab3b"   "Rgs3"    "Rnase4"  "Rprm"    "Sfrp1"   "Socs1"   "Sprr2a3"
  # [36] "Susd4"   "Tapbpl"  "Tgfb2"   "Zbp1"  
  (n.GS <- length(GS.CTS))  # 39
  write.table(GS.CTS, file=paste0('GS.CTS.CMP_at.least.',n.select,'identificaitons_fr.',length(select),'predicitons.txt'),
              row.names=F, col.names=F)
  
  ## 2.5) evaluate the F1 score for the identified CTS genes 
  ## at mesoderm-endoderm bifurcation and MP, respectively   
  ## using the proxy gold standard ----------------------------------
  # load the CTS prediction with different clustering methods
  running.methods <- list.dirs()
  running.methods <- running.methods[grepl("./C_", running.methods)]
  running.methods
  methods <- lapply(running.methods, function(x) unlist(strsplit(x, "./C_"))[2]) %>% unlist()
  x <- which(is.na(methods))
  if(length(x)>0) {
    methods <- methods[-x]
    running.methods <- running.methods[-x]
  }
  ## manually adding the QuanTC's predictions
  running.methods <- c(running.methods, './QuanTC_run')
  methods <- c(methods,'QuanTC_run')
  
  which.CT <- c('MEP','MP','CMP')
  fid <- c('GS.CTS.MEP_at.least.3identificaitons_fr.5predicitons.txt',
           'GS.CTS.MP_at.least.4identificaitons_fr.9predicitons.txt',
           'GS.CTS.CMP_at.least.3identificaitons_fr.6predicitons.txt')  
  
  # F1.res <- get.multiple.F1.score(fid, match.col.name=which.CT, 
  #                                 QuanTC.genes.list.id = 1:6,
  #                                 running.methods=running.methods, methods=methods,
  #                                 sce=sce,
  #                                 best.matched.clusterID=best.matched.clusterID)
  F1.res <- get.multiple.F1.score(match.col.name='MEP', 
                                  CTS.reference = read.table(fid[1])[,1],
                                  QuanTC.genes.list.id = 1:6,
                                  running.methods=running.methods, methods=methods,
                                  sce = sce,
                                  best.matched.clusterID=best.matched.clusterID)
  F1.KLF1 <- F1.res$F1.scores$MEP
  F1.KLF1.ctl <- F1.res$F1.ctl$MEP
  round(F1.res$F1.scores$MEP,3)
  # consensus_ks13 consensus_ks15 consensus_ks17 consensus_ks19     Leiden_0.4 
  #    0.077          0.238          0.000          0.000          0.000 
  # Leiden_0.8     Leiden_1.2   SNNGraph_k10    SNNGraph_k5           Soft 
  #    0.480          0.000          0.516          0.485          0.000 
  # Soft_wo.TC     QuanTC_run 
  #    0.000         0.028
  
  F1.res <- get.multiple.F1.score(match.col.name='MP', 
                                  CTS.reference = read.table(fid[2])[,1],
                                  QuanTC.genes.list.id = 1:6,
                                  running.methods=running.methods, methods=methods,
                                  sce = sce,
                                  best.matched.clusterID=best.matched.clusterID)
  F1.MP <- F1.res$F1.scores$MP
  F1.MP.ctl <- F1.res$F1.ctl$MP
  round(F1.res$F1.scores$MP,3)
   # consensus_ks13 consensus_ks15 consensus_ks17 consensus_ks19     Leiden_0.4 
   # 0.290          0.347          0.238          0.250          0.394 
   # Leiden_0.8     Leiden_1.2   SNNGraph_k10    SNNGraph_k5           Soft 
   # 0.453          0.430          0.073          0.163          0.000 
   # Soft_wo.TC     QuanTC_run 
   # 0.000          0.000 
  
  F1.res <- get.multiple.F1.score(match.col.name='CMP', 
                                  CTS.reference = read.table(fid[3])[,1],
                                  QuanTC.genes.list.id = 1:6,
                                  running.methods=running.methods, methods=methods,
                                  sce = sce,
                                  best.matched.clusterID=best.matched.clusterID)
  F1.CMP <- F1.res$F1.scores$CMP
  F1.CMP.ctl <- F1.res$F1.ctl$CMP
  round(F1.res$F1.scores$CMP,3)
  # consensus_ks13 consensus_ks15 consensus_ks17 consensus_ks19     Leiden_0.4 
  # 0.246          0.297          0.000          0.000          0.000 
  # Leiden_0.8     Leiden_1.2   SNNGraph_k10    SNNGraph_k5           Soft 
  # 0.117          0.280          0.000          0.418          0.116 
  # Soft_wo.TC     QuanTC_run 
  # 0.238          0.000  
   
  ## calculating F1 scores for the gene module with the highest DNB score
  # which refer to the SNNGraph_k5 run, but the C11 did not have the highest DNB score 
 
  nn = length(F1.KLF1)  # 12
  Normalized.F1 <- Normalize.F1(c(F1.KLF1.ctl, F1.KLF1, F1.MP.ctl, F1.MP, F1.CMP.ctl, F1.CMP))
  Norm.F1.KLF1.ctl <- Normalized.F1[1:nn] 
  Norm.F1.KLF1 <- Normalized.F1[(nn+1):(2*nn)] 
  Norm.F1.MP.ctl <- Normalized.F1[(2*nn+1):(3*nn)]
  Norm.F1.MP <- Normalized.F1[(3*nn+1):(4*nn)]  
  Norm.F1.CMP.ctl <- Normalized.F1[(4*nn+1):(5*nn)]
  Norm.F1.CMP <- Normalized.F1[(5*nn+1):(6*nn)]  
  
}

## Verify if the CTS identified the gene of interest which is Etv2                                               ##
## only the 'SNNGraph_allcells' did                                                                                                              ##
###################################################################################################################
# load the CTS prediction with different clustering methods
HEP.marker <-  c('Etv2','Rhoj','Tal1','Lyl1','Plvap','Rasip1')   

HEP.CTS <- NULL 
for(i in 1:length(running.methods)){
  if(grepl('QuanTC', methods[i] )) {
    load(file=paste0(running.methods[i],'/CTS.RData'))
    set2 <- CTS 
    tmp <- HEP.marker %in% unlist(set2) 
    HEP.CTS = rbind(HEP.CTS, tmp)
  } else {
    load(file=paste0(running.methods[i],'/BioTIP.res.RData'))
    if(length(res)>1){ # flag 1
      names(res)
      #[1] "CTS.candidate" "CTS.score"     "Ic.shrink"     "significant" 
      set2 <- res$CTS.candidate[which(res$significant)]
      
      best.match <- best.matched.clusterID[paste0('C_',methods[i]),'MEP']
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
HEP.CTS
#                 Etv2    Rhoj    Tal1    Lyl1    Plvap   Rasip1 
# consensus_ks13 "FALSE" "FALSE" "FALSE" "FALSE" "FALSE" "FALSE"
# consensus_ks15 "FALSE" "TRUE"  "TRUE"  "FALSE" "FALSE" "FALSE"
# consensus_ks17 "FALSE" "FALSE" "FALSE" "FALSE" "FALSE" "FALSE"
# consensus_ks19 "FALSE" "FALSE" "FALSE" "FALSE" "FALSE" "FALSE"
# Leiden_0.4     "FALSE" "FALSE" "FALSE" "FALSE" "FALSE" "FALSE"
# Leiden_0.8     "TRUE"  "TRUE"  "TRUE"  "FALSE" "FALSE" "FALSE"
# Leiden_1.2     "FALSE" "TRUE"  "FALSE" "FALSE" "FALSE" "TRUE" 
# SNNGraph_k10   "TRUE"  "TRUE"  "TRUE"  "FALSE" "FALSE" "FALSE"
# SNNGraph_k5    "TRUE"  "TRUE"  "TRUE"  "FALSE" "FALSE" "FALSE"
# Soft           "FALSE" "FALSE" "FALSE" "FALSE" "FALSE" "FALSE"
# Soft_wo.TC     "FALSE" "FALSE" "FALSE" "FALSE" "FALSE" "FALSE"
# QuanTC_run     "FALSE" "FALSE" "TRUE"  "FALSE" "FALSE" "FALSE"


###########################################################
## Section 3) reproting table, compare BioTIP to QuanTC  ##
###########################################################
{
  tb <- data.frame(Jaccard.CT = round(jaccard.CT,3), 
                   F1.KLF1.ctl = round(F1.KLF1.ctl,3), 
                   F1.KLF1 = round(F1.KLF1,3),
                   F1.MP.ctl = round(F1.MP.ctl,3), 
                   F1.MP = round(F1.MP,3),
                   F1.CMP.ctl = round(F1.CMP.ctl,3), 
                   F1.CMP = round(F1.CMP,3),
                   Norm.F1.KLF1.ctl = round(Norm.F1.KLF1.ctl,3), 
                   Norm.F1.KLF1 = round(Norm.F1.KLF1,3),
                   Norm.F1.MP.ctl = round(Norm.F1.MP.ctl,3), 
                   Norm.F1.MP = round(Norm.F1.MP,3),
                   Norm.F1.CMP.ctl = round(Norm.F1.CMP.ctl,3), 
                   Norm.F1.CMP = round(Norm.F1.CMP,3)
                   )
  nrow(tb) #[1] 12
  
  # reorder to show
  myorder <- c('C_SNNGraph_k5','C_SNNGraph_k10', 
               paste0('C_consensus_ks',c(13,15,17)),
               paste0('C_Leiden_',c(0.4, 0.8, 1.2)),
               "C_Soft_wo.TC","C_Soft",'QuanTC_run' 
  )
  tb$detected.CT <- lapply(CT, toString) %>% unlist()
  
  # (mytable <- tb[myorder,
  #                c('detected.CT','Jaccard.CT', 'F1.MEP', 'F1.nor.MEP','F1.MP', 'F1.nor.MP',
  #                  'F1.CMP', 'F1.nor.CMP')])
  #                                    detected.CT Jaccard.CT F1.MEP F1.nor.MEP F1.MP F1.nor.MP F1.CMP F1.nor.CMP
  # C_SNNGraph_k5             5, 3, 6, 10, 11, 4, 7      0.429  0.485      0.937 0.163     0.370  0.418      0.968
  # C_SNNGraph_k10                 9, 6, 2, 1, 4, 3      0.286  0.516      0.953 0.073     0.195  0.000      0.168
  # C_consensus_ks13      10, 5, 12, 6, 4, 7, 13, 9      0.375  0.077      0.370 0.290     0.660  0.246      0.758
  # C_consensus_ks15 7, 14, 12, 4, 15, 11, 6, 13, 8      0.333  0.238      0.657 0.347     0.772  0.297      0.851
  # C_consensus_ks17         10, 3, 7, 14, 8, 16, 2      0.111  0.000      0.247 0.238     0.542  0.000      0.168
  # C_Leiden_0.4                            6, 5, 4      0.200  0.000      0.247 0.394     0.846  0.000      0.168
  # C_Leiden_0.8           12, 4, 10, 9, 8, 3, 6, 1      0.375  0.480      0.935 0.453     0.914  0.117      0.431
  # C_Leiden_1.2          13, 10, 3, 11, 8, 6, 9, 7      0.222  0.000      0.247 0.430     0.891  0.280      0.823
  # C_Soft_wo.TC                     C7, C8, C5, C6      0.167  0.000      0.247 0.000     0.099  0.238      0.741
  # C_Soft                                   C8, C5      0.250  0.000      0.247 0.000     0.099  0.116      0.428
  # QuanTC_run       C1, C4, C2, C6, C5, C8, C3, C7      0.222  0.000      0.247 0.000     0.099  0.000      0.168
  
  write.table(tb, 
              file='jaccard.CT_F1.CTS.txt', 
              sep='\t')
  
  
  ###############
  ## PLOT #######
  ###############
  
  ## generate bar plot #######################
  df <- read.table('jaccard.CT_F1.CTS.txt', header=T, sep='\t')
  tmp <- rownames(df)
  tmp = gsub('C_', '', tmp)
  #df$method = gsub('QuanTBioTIP_run', 'QuanTC_Soft', tmp)
  df$method=tmp
  df$method=factor(df$method, levels=rev(df$method))
  
  p1 <- ggbarplot(df, "method", "Jaccard.CT", orientation = "horiz") +
    geom_col(aes(fill = Jaccard.CT)) + 
    scale_fill_gradient2(low = "white",high = "green") +
    theme_bw(base_size=10)  +        # use a black-and-white theme with set font size
    theme(legend.position = "none") + coord_flip()   
  p2 <- ggbarplot(df, "method", "Norm.F1.KLF1", orientation = "horiz") +
    geom_col(aes(fill = Norm.F1.KLF1)) + 
    scale_fill_gradient2(low = "white",high = "blue") +
    theme_bw(base_size=10)    +      # use a black-and-white theme with set font size
    theme(legend.position = "none")  + coord_flip()  
  
  p3 <- ggbarplot(df, "method", "Norm.F1.MP", orientation = "horiz") +
    geom_col(aes(fill = Norm.F1.MP)) + 
    scale_fill_gradient2(low = "white",high = "blue") +
    theme_bw(base_size=10)    +      # use a black-and-white theme with set font size
    theme(legend.position = "none")  + coord_flip()  
  p4 <- ggbarplot(df, "method", "Norm.F1.CMP", orientation = "horiz") +
    geom_col(aes(fill = Norm.F1.CMP)) + 
    scale_fill_gradient2(low = "white",high = "blue") +
    theme_bw(base_size=10)    +      # use a black-and-white theme with set font size
    theme(legend.position = "none")  + coord_flip()  
  
  pdf(file='bar_robustness.pdf', width=12, height=4)
  gridExtra::grid.arrange(
    p1, p2, p3, p4 
    ,ncol=4)
  dev.off()
  }
  
  
###############################################################
## Section 4) stability using the SNNGraph (k=5) clusters    ##
###############################################################
load(file="ShrinkM.RData")
dim(M) #4000 4000

colnames(colData(sce))
# [1] "Sample"           "Barcode"          "sizeFactor"       "label"           
# [5] "C_SNNGraph_k5"    "C_SNNGraph_k10"   "C_consensus_ks13" "C_consensus_ks15"
# [9] "C_consensus_ks17" "C_consensus_ks19" "C_Leiden_0.4"     "C_Leiden_0.8"    
# [13] "C_Leiden_1.2"     "C_Soft"           "QuanTC_run"  
sce$label <- sce$C_SNNGraph_k5  # selected clustering method !!!
rowData(sce) <- NULL # to allow run sample()

######### 4.1) setting parameters, the same as section 1  ---------------------------
## but try multiple getTopMCI.gene.minsize 
getTopMCI.gene.minsize <- c(20,30,40)

######### 4.2) run BioTIP on downsized data --------------------
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
                              getTopMCI.gene.minsize=getTopMCI.gene.minsize[m], 
                              MCIbottom=MCIbottom[i],
                              localHVG.preselect.cut=localHVG.preselect.cut,  #logmat.local.HVG.testres=testres,
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
    save(res, flag, counts, file=paste0('stability/BioTIP_res_20runs_MCIbottom',MCIbottom[i],
                                        '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    
  }
}



## 4.3) evaluate the stability of CT detection
## using the robustly detected C2 as reference -----------------------------------
getTopMCI.gene.minsize
#[1] 20 30 40
MCIbottom 
#[1] 1.5 2.0
{
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
  save(CT.list, file= paste0('stability/CT.list_variable.MCIbottom_Modulesize',
                             getTopMCI.gene.minsize[m],'.Rdata'))
  
  df <- data.frame(Frequency = as.numeric(Freq),
                   clusterID = rep(rownames(Freq), ncol(Freq)),
                   MCIbottom = lapply(as.character(MCIbottom ), function(x) rep(x,nrow(Freq))) %>% unlist()
  )
  head(df)
  df$clusterID <- factor(df$clusterID, levels =levels(sce$label)) #!!!
  df$moduleSize <- rep(getTopMCI.gene.minsize[m], nrow(df))
  
  df.CT <- rbind(df.CT, df)  
}
df.CT$Frequency = df.CT$Frequency/n.scability

save(df.CT, 
     file= 'stability/CT.detection_variable.MCIbottom_Modulesize.Rdata')
}
head(df.CT)
#     Frequency clusterID MCIbottom moduleSize
# 1      0.10         1       1.5         20
# 2      1.00        10       1.5         20
# 3      0.60        11       1.5         20
# 4      0.70        12       1.5         20
# 5      0.85        14       1.5         20
# 6      0.00        17       1.5         20


## 4.3) plot for variable Modulesize, when MCIbottom==2         ##
#############################################################
load(file= 'stability/CT.detection_variable.MCIbottom_Modulesize.Rdata')
df.CT$moduleSize <- as.character(df.CT$moduleSize )

ggplot(data=subset(df.CT, MCIbottom==2),
       aes(x=clusterID, y=Frequency, group=moduleSize)) +
  geom_line(aes(color=moduleSize), size=1.5)+
  geom_point(aes(color=moduleSize)) 
dev.copy2pdf(file='stability/CT.detection_MCIbottom1.5_variable.Modulesize.pdf')

ggplot(data=subset(df.CT, MCIbottom==1.5),
       aes(x=clusterID, y=Frequency, group=moduleSize)) +
  geom_line(aes(color=moduleSize), size=1.5)+
  geom_point(aes(color=moduleSize)) 
dev.copy2pdf(file='stability/CT.detection_MCIbottom2_variable.Modulesize.pdf')


## 4.4) evaluate the stability of CTS identification
## for the MBP and MP, respectively
###########################################################
load(file= paste0('stability/CT.list_variable.MCIbottom_Modulesize',
                  getTopMCI.gene.minsize[2],'.Rdata'))
names(CT.list)
#[1] "1.5" "2" 


# load('original_run/CTS.RData')
# names(CTS)
# #[1] "9"  "10" "9" 
# C9.ref = CTS[which(names(CTS)=='9')]
# C10.ref = CTS['10']
MEP.ref <- read.table('GS.CTS.MEP_at.least.3identificaitons_fr.5predicitons.txt')
MP.ref <- read.table( 'GS.CTS.MP_at.least.4identificaitons_fr.9predicitons.txt')  
CMP.ref <- read.table( 'GS.CTS.CMP_at.least.3identificaitons_fr.6predicitons.txt')  
ref <- list('11'= MEP.ref[,1],'6'= MP.ref[,1], '4'= CMP.ref[,1])

MCIbottom  # 1.5,  2
i=2
for(k in 1:length(ref)){
  CT.ref <- ref[k]
  
  F1.CT <- F1.CT.ctl <- N.module <- list()
  for(m in 1:length(getTopMCI.gene.minsize)){
    load(file=paste0('stability/BioTIP_res_20runs_MCIbottom',MCIbottom[i],
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
       file=paste0('stability/F1.CTS_MCIbottom2_GS.C',names(ref)[k],'.RData'))
}

load('stability/F1.CTS_MCIbottom2_GS.C11.RData')
F1.MEP <- F1.CT 
F1.MEP.ctl <- F1.CT.ctl
lapply(N.module, function(x) mean(x, na.rm=T)) %>% unlist()
#      20       30       40 
# 58.54545 65.52941 63.26667  
load('stability/F1.CTS_MCIbottom2_GS.C6.RData')
F1.MP <- F1.CT 
F1.MP.ctl <- F1.CT.ctl
lapply(N.module, function(x) mean(x, na.rm=T)) %>% unlist()
#       20       30       40 
# 75.10000 75.25000 73.21053 
load('stability/F1.CTS_MCIbottom2_GS.C4.RData')
F1.CMP <- F1.CT 
F1.CMP.ctl <- F1.CT.ctl
lapply(N.module, function(x) mean(x, na.rm=T)) %>% unlist()
#       20       30       40 
#35.90000 46.66667 57.50000 


lapply(F1.MEP, mean) %>% unlist()
#       20        30        40
# 0.1909455 0.3166479 0.2853768 
lapply(F1.MP, mean) %>% unlist()
#       20        30        40 
# 0.2052842 0.2040435 0.1786371 

df <- data.frame(F1 = c(unlist(F1.MEP.ctl),unlist(F1.MEP),
                        unlist(F1.MP.ctl),unlist(F1.MP),
                         unlist(F1.CMP.ctl),unlist(F1.CMP) ),
                 minModuleSize = c(rep(getTopMCI.gene.minsize, lengths(F1.MEP.ctl)),
                                   rep(getTopMCI.gene.minsize, lengths(F1.MEP)),
                                   rep(getTopMCI.gene.minsize, lengths(F1.MP.ctl)),
                                   rep(getTopMCI.gene.minsize, lengths(F1.MP)),
                                   rep(getTopMCI.gene.minsize, lengths(F1.CMP.ctl)),
                                   rep(getTopMCI.gene.minsize, lengths(F1.CMP)) ),
                 CT = c(rep('MEP', sum(lengths(F1.MEP.ctl))+sum(lengths(F1.MEP))),
                        rep('MP', sum(lengths(F1.MP.ctl))+sum(lengths(F1.MP))),
                        rep('CMP', sum(lengths(F1.CMP.ctl))+sum(lengths(F1.CMP))) ),
                 CTS = c(rep('ctl', sum(lengths(F1.MEP.ctl))), rep('CTS', sum(lengths(F1.MEP))),
                         rep('ctl', sum(lengths(F1.MP.ctl))), rep('CTS', sum(lengths(F1.MP))),
                         rep('ctl', sum(lengths(F1.CMP.ctl))), rep('CTS', sum(lengths(F1.CMP))) ) 
)
df$minModuleSize <- as.character(df$minModuleSize) 
df$CTS <- factor(df$CTS, levels = c('ctl','CTS')) 
df$color.lab <- paste0(df$minModuleSize,df$CTS)
df$color.lab <- factor(df$color.lab, levels=c('20ctl','20CTS','30ctl','30CTS','40ctl','40CTS')) #, '40ctl','40CTS'))
df$CT <- factor(df$CT, levels =c('MEP','MP','CMP'))

df$Norm.F1 = Normalize.F1(df$F1)


library(ggbeeswarm)
library(grid)
library(gridExtra)
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

pdf(file='stability/CTS.F1_variable.minModuleSize.pdf')
gridExtra::grid.arrange(
  p1, p2, top = "CTS stability", bottom = "different module size settings"
  ,ncol=1)
dev.off()

dim(sce)  # 4000 1531
## prsent the statistics using proxy gold standard here an in the figure
tmp <- compare_means(Norm.F1 ~ color.lab, data=df, group='CT', method='t.test')
subset(tmp, CT=='MEP' & group1=='20ctl' & group2=='20CTS')
#1 MEP   Norm.F1 20ctl  20CTS  0.000411 0.0082 0.00041  ***      T-test
subset(tmp, CT=='MEP' & group1=='30ctl' & group2=='30CTS')
# 1 MEP   Norm.F1 30ctl  30CTS  0.00000000407 0.00000012 4.1e-09  ****     T-test
subset(tmp, CT=='MEP' & group1=='40ctl' & group2=='40CTS')
#1 MEP   Norm.F1 40ctl  40CTS  0.000000499 0.000013 5.0e-07  ****     T-test

tmp <- compare_means(Norm.F1 ~ color.lab, data=df, group='CT', method='t.test')
subset(tmp, CT=='MP' & group1=='20ctl' & group2=='20CTS')
# MP    Norm.F1 20ctl  20CTS  4.61e-14 2e-12 4.6e-14  ****     T-test
subset(tmp, CT=='MP' & group1=='30ctl' & group2=='30CTS')
# 1 MP    Norm.F1 30ctl  30CTS  9.09e-13 3.7e-11 9.1e-13  ****     T-test
subset(tmp, CT=='MP' & group1=='40ctl' & group2=='40CTS')
#1 MP    Norm.F1 40ctl  40CTS  1.83e-10 0.0000000066 1.8e-10  ****     T-test

tmp <- compare_means(Norm.F1 ~ color.lab, data=df, group='CT', method='t.test')
subset(tmp, CT=='CMP' & group1=='20ctl' & group2=='20CTS')
#1 CMP   Norm.F1 20ctl  20CTS  1.75e-10 0.0000000065 1.7e-10  ****     T-test
subset(tmp, CT=='CMP' & group1=='30ctl' & group2=='30CTS')
#1 CMP   Norm.F1 30ctl  30CTS  0.00000000285 0.000000091 2.9e-09  ****     T-test
subset(tmp, CT=='CMP' & group1=='40ctl' & group2=='40CTS')
#1 CMP   Norm.F1 40ctl  40CTS  0.00000782 0.00019 7.8e-06  ****     T-test


