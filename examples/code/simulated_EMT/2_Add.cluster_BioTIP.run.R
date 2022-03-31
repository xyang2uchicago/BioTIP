# this code does the following analyses
# 1) load downloaded data matrix to an R object and cluster the cells using different clustering methods
# 2) run BioTIP analysis on the data with different cell clusterings
# 3) assess the jacaja similarity of the predicted CT clusters, the F1 score for the predicted CTSs   
# 4) 

setwd('F:/projects/BioTIP/result/QuanTC_simulation/')


# the following funciton has been updated on github, but before updating the package to Bioconductor, you need to call it locally to test
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/ForJennifer/optimize.sd_selection.R') 
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/BioTIP.wrap.R')
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/stability_score.R') 

library(BioTIP)
library(RColorBrewer)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(scater)
library(scran)
library(pheatmap)

library(dplyr)

library(ggplot2)
library(gridExtra)  # required by grid.arrange()
library(ggpubr)  # required by p-values on bar

require(stringr)
require(psych)
require(igraph)

library(SC3)
packageVersion('SC3') # 1.18.0


######################################################################
## Section 1)                                                       ##
## translate the MATHLAB data into an R object                      ##
## original code: 2_QuanTC_simulation_cluster_notused               ##
## First cluster all cells using different methods and parameters   ##
## Second cluster not TC cells only                                 ##
## DO NOT repeat                                                    ##
######################################################################
{
  #################################################################
  ####### load profiles ######################
  # refer to 1_export_data.m
  
  X <- NULL # a matrix of 18 gene x 5363 cells
  for(i in 1:5){
    fileName = paste0(inputpath,'data',i,'.txt')
    tmp <- read.table(fileName, sep=",")
    dim(tmp) # 312 18,  a cell x gene matrix
    X <- cbind(X, t(tmp))
  } 
  
  # check normalitity
  # boxplot(X)
  # boxplot(log2(X+1))
  hist(log2(X+1), 50, ylab='log2(x+1)', main="data of 5363 cells")
  dev.copy2pdf(file='hist_data_18gene_exp.pdf')
  
  y <- NULL # a matrix of 18 gene x 5363 cells
  for(i in 1:5){
    fileName = paste0(inputpath,'true_label',i,'.txt')
    tmp <- read.table(fileName, sep=",")
    dim(tmp) # 312 1 
    y <- c(y, tmp[,1])
  } 
  length(y)
  # [1] 5363
  table(y)
  #   1    2    3    4    5 
  # 577 1849  801  985 1151 
  tmp = y
  tmp[tmp==1] = 'E'
  tmp[tmp==2] = 'I1'
  tmp[tmp==3] = 'I2'
  tmp[tmp==4] = 'M'
  tmp[tmp==5] = 'TC'
  tmp = factor(tmp, levels = c('E','I1', 'I2', 'M', 'TC'))
  y = factor(y, levels=c('1','2','3','4','5'))
  
  fileName = paste0(inputpath,'gene_name.txt')
  gene_name <- read.table(fileName, sep=",")
  gene_name=as.character(gene_name[1,])
  
  rownames(X) <- gene_name
  colnames(X) <- y
  
  sce <- SingleCellExperiment(list(counts=X, logcounts=log2(X+1)),
                              colData=DataFrame(true_label=y,
                                                cell_type=tmp),
                              rowData=DataFrame(gene_name=gene_name),
                              metadata=list(study="EMT simulation genearted by QuanTC authors")
  )
  sce.all <- sce
  
  DBs <- c('all','woTC')
  #for(i in 1:length(DBs)){
  m=1 {
    if(DBs[i]=='all') sce <- sce.all else sce <- sce.all[,-which(sce$cell_type=='TC')]
   
  #################################################################
  ## reduce dimention 
  
  library(scater)
  set.seed(1111001)
  sce <- scater::runPCA(sce,  exprs_values='logcounts') 
  
  # choose the top 10 PCs when working on log2counts
  percent.var <- attr(reducedDim(sce), "percentVar")
  plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
  fileName <- ifelse(methods[i]=='all','PC_varnace.pdf','PC_varnace.wo.TC.pdf' )
  dev.copy2pdf(file=fileName)
  
  
  library(scater)
  set.seed(1111001)
  sce <- runTSNE(sce, dimred='PCA')
  set.seed(1111001)
  sce <- runUMAP(sce, dimred='PCA')
  
  reducedDimNames(sce)
  #[1] "PCA"           "TSNE"         "UMAP"  
  
  
  ## remove the repeated sample names, just need to run once
  if(i==1) colnames(sce) = 1:dim(sce)[2]
  ##########################################
  ## 1.1) SNNGraph clustering method  ######
  ##########################################
  library(scran)
  g.200 <- buildSNNGraph(sce, k=200, use.dimred = 'PCA')
  clust <- igraph::cluster_walktrap(g.200)$membership
  table(clust)
  #   1    2    3    4 
  # 2373 1017  752 1221
  sce$C_SNNGraph.k200 = factor(clust)
  
  g.100 <- buildSNNGraph(sce, k=100, use.dimred = 'PCA')
  clust <- igraph::cluster_walktrap(g.100)$membership
  table(clust)
  #   1    2    3    4    5 
  # 1798  885  711  753 1216
  sce$C_SNNGraph.k100 = factor(clust)
  
  g.20 <- buildSNNGraph(sce, k=20, use.dimred = 'PCA')
  clust <- igraph::cluster_walktrap(g.20)$membership
  table(clust)
  #   1    2    3    4    5    6    7    8
  # 536  837 1238  782  539  931  214  286
  sce$C_SNNGraph.k20 = factor(clust)
  
  

 
  ################################################################
  # 13) # extract Leiden clustering (using Seurat)  
  # https://satijalab.org/seurat/articles/get_started.html
  # Leiden requires the leidenalg python.
  # We apply the Seurat packge, by setting the algorithm =4 in the function FindClusters()	
  # This parameter decides the algorithm for modularity optimization 
  # (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 
  # 3 = SLM algorithm; 4 = Leiden algorithm). 
  # The resolution parameter: use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
  # Seurat author recommended that:
  # We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
  ###################################################################
   # convert from SingleCellExperiment
  sce.seurat <- as.Seurat(sce)
  # Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
  sce.seurat
  # An object of class Seurat 
  # 96 features across 929 samples within 1 assay 
  # Active assay: originalexp (96 features, 0 variable features)
  # 3 dimensional reductions calculated: PCA, TSNE, UMAP
  
  
  # Computes the k.param nearest neighbors for a given dataset
  sce.seurat <- FindNeighbors(sce.seurat, reduction = "PCA", k.param = 20) # k.param = 20 is default
  
  sce.seurat <- FindClusters(sce.seurat, resolution = 0.4, algorithm = 4) # smaller  number of communities 
  table(Idents(sce.seurat), as.vector(colData(sce)$cell_type))  
  #     E   I1   I2    M   TC
  # 1    0 1125    2    0    9
  # 2    0    0    0  946   20
  # 3    0   25  753    0   21
  # 4    0   44   38   28  683
  # 5  577    0    0    0  173
  # 6    0  655    0    0   14
  # 7    0    0    8   11  231
  colData(sce)$C_Leiden_0.4 = Idents(sce.seurat)
  
  sce.seurat <- FindClusters(sce.seurat, resolution = 0.8, algorithm = 4) # smaller  number of communities 
  table(Idents(sce.seurat), as.vector(colData(sce)$cell_type))  
  #      E  I1  I2   M  TC
  # 1    0   0   0 919  17
  # 2    0  44  37   9 625
  # 3    0 558   2   0   6
  # 4    0 439   0   0   0
  # 5    0   1 419   0  10
  # 6    0 422   0   0   7
  # 7    0  24 335   0  13
  # 8    0 361   0   0   9
  # 9    0   0   8  38 235
  # 10 270   0   0   0   7
  # 11 271   0   0   0   4
  # 12  36   0   0   0 162
  # 13   0   0   0  19  56
  
  colData(sce)$C_Leiden_0.8 = Idents(sce.seurat)
  
  sce.seurat <- FindClusters(sce.seurat, resolution = 1.2, algorithm = 4) # smaller  number of communities 
  table(Idents(sce.seurat), as.vector(colData(sce)$cell_type))  
  #      E  I1  I2   M  TC
  # 1    0   0   0 520  12
  # 2    0 462   0   0   4
  # 3    0 452   2   0   5
  # 4    0 444   0   0   0
  # 5    0   0   0 404   5
  # 6    0  24 332   0  12
  # 7    0  18   5   0 332
  # 8  270   0   0   0   7
  # 9    0   0   8  33 235
  # 10 271   0   0   0   4
  # 11   0 271   0   0   2
  # 12   0   1 262   0   8
  # 13   0  23   2   0 206
  # 14   0   1  31  28 146
  # 15  36   0   0   0 162
  # 16   0 153   0   0  10
  # 17   0   0 159   0   1
  colData(sce)$C_Leiden_1.2 = Idents(sce.seurat)
  
 
  sce
  # class: SingleCellExperiment 
  # dim: 18 5363 
  # metadata(1): study
  # assays(2): counts logcounts
  # rownames(18): snailt SNAIL ... Ncad Ovol2
  # rowData names(1): gene_name
  # colnames(5363): 1 2 ... 5362 5363
  # colData names(8): true_label cell_type ... C_Leiden_0.8 C_Leiden_1.2
  # reducedDimNames(3): PCA TSNE UMAP
  # altExpNames(0):
  
  if(DBs[i]=='all') save(sce, file='sce.RData') else save(sce, file='sce.wo.TC.Rdata')
  
  }
}


##########################################
## section 2)                           ##
##    BioTIP analysis                   ##
## for all cells                        ##
##########################################
DBs #[1] "all"  "woTC"


#for(i in 1:length(DBs)){
m=1 {
  if(DBs[i]=='all') setwd('F:/projects/BioTIP/result/QuanTC_simulation/simulatedEMT_robustness/') else
    setwd('F:/projects/BioTIP/result/QuanTC_simulation/simulatedEMT.wo.TC_robustness/')
  
  if(DBs[i]=='all') load('../sce.RData') else load('../sce.wo.TC.RData')
  
  ######### 1) setting parameters,   ######################
  local.HVG.optimize = FALSE # shutoff the optimazation process for only 18 genes
  localHVG.preselect.cut = 1  # use all 18 genes
  getNetwork.cut.fdr = 0.05  # to construct RW network and extract co-expressed gene moduels
  getTopMCI.gene.minsize = 6  # min number of genes in an identified CTS
  getTopMCI.n.states = 2  # A number setting the number of states to check the significance of DNB score (i.e. MCI) 
  # This parameter can be enlarge when inputting more clusters
  MCIbottom = 0  # A number setting the threshold to prioritize the initially selected CTS candidates.
  # In our experiment, a number between 2 and 4 generated expected resutls for real scRNA-seq datasets.
  # set to 0 here because there are only 18 genes tested !!
  
  # # M is precalculated correlation matrix, will be reused in the downstream simulation analysis
  logmat <- as.matrix(logcounts(sce))
  M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
  dim(M) # 18  18
  save(M, file="CTS_ShrinkM.RData", compress=TRUE) # save the file as its calculation runs a while
  rm(logmat)
  
  ############### all new clustering outputs  ##################
  x <- grep("C_", colnames(colData(sce)))
  colnames(colData(sce))[x]
  # [1] "C_SNNGraph.k200" "C_SNNGraph.k100" "C_SNNGraph.k20"  "C_Leiden_0.4"    "C_Leiden_0.8"   
  # [6] "C_Leiden_1.2"
  
  # set.seed(102020)
  for(i in 1:length(x)){
    samplesL <- split(rownames(colData(sce)), f = colData(sce)[x[i]])
    (tmp = lengths(samplesL))
    
    res <- BioTIP.wrap(sce, samplesL, subDir=colnames(colData(sce))[x[i]],
                       getTopMCI.n.states=getTopMCI.n.states, 
                       getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                       MCIbottom=MCIbottom,
                       localHVG.preselect.cut=localHVG.preselect.cut, 
                       getNetwork.cut.fdr=getNetwork.cut.fdr, 
                       M=M,
                       empirical.MCI.p.cut=1.1, ## because there are only 18 genes tested !!
                       empirical.IC.p.cut = 1.1, ## because there are only 18 genes tested !!
                       permutation.method='both',
                       verbose=TRUE, plot=TRUE)           
    save(res, file=paste0(colnames(colData(sce))[x[i]],'/BioTIP.res.RData'))  
    
  }
  
  ## test for the true label given by the author (cell_type)
  samplesL <- split(rownames(colData(sce)), f = sce$cell_type)
  (tmp = lengths(samplesL))
  
  res <- BioTIP.wrap(sce, samplesL, subDir='C_cell_type',
                     getTopMCI.n.states=getTopMCI.n.states, 
                     getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                     MCIbottom=MCIbottom,
                     localHVG.preselect.cut=localHVG.preselect.cut, 
                     getNetwork.cut.fdr=getNetwork.cut.fdr, 
                     M=M,
                     empirical.MCI.p.cut=1.1, ## because there are only 18 genes tested !!
                     empirical.IC.p.cut = 1.1, ## because there are only 18 genes tested !!
                     permutation.method='both',
                     verbose=TRUE, plot=TRUE)           
  save(res, file='C_cell_type/BioTIP.res.RData')  
  
}



################################################################
## section 3)                                                 ##
## assess the jacaja similarity of the predicted CT clusters  ##
## assess the F1 score for the predicted CTSs                 ##
################################################################

# source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/stability_score.R') 
# source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/BioTIP.wrap.R')
# 
# 
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


#for(m in 1:length(DBs))
m=1 {
  if(DBs[m]=='all') setwd('F:/projects/BioTIP/result/QuanTC_simulation/simulatedEMT_robustness/') else
    setwd('F:/projects/BioTIP/result/QuanTC_simulation/simulatedEMT.wo.TC_robustness/')
  
  if(DBs[m]=='all') load('../sce.RData') else load('../sce.wo.TC.RData')
  
  ## 3.1) Summarize the CT identifications for the CT cells of true label  ##########
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
  
  ## manually add the QuanTC's prediction about CT states is M->I1
  ## See F:\projects\BioTIP\result\QuanTC_simulation\QuanTC for details
  CT$QuanTC_run <- c('TC','I1')
  CT
  
 ## 3.2) Annotate cell identities by comparing to the original cluster IDs. ##
  colnames(colData(sce))
  # [1] "true_label"      "cell_type"       "C_SNNGraph.k200" "C_SNNGraph.k100"
  # [5] "C_SNNGraph.k20"  "C_Leiden_0.4"    "C_Leiden_0.8"    "C_Leiden_1.2"   
  # [9] "C_SNNGraph.k50" 
  x <- grep("C_", colnames(colData(sce)))
  colnames(colData(sce))[x]
  
  best.matched.clusterID <- data.frame()
  query = c('TC','I1') # the cluster ID for bifurcation at PS, CM, EP, and Ecoderm
  for(i in x){
    tmp <- table(sce$cell_type, colData(sce)[,i])
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
  colnames(best.matched.clusterID) <-  query
  
  ## add QuanTC's identification
  best.matched.clusterID <- rbind(best.matched.clusterID, 
                                  'QuanTC_run'=  query)
  
  setdiff(names(CT), rownames(best.matched.clusterID))
  # [1] "C_cell_type"
  
  # Now, manually add the best-match infor for these additional BioTIP runs
  best.matched.clusterID <- rbind(best.matched.clusterID, 
                                  'C_cell_type'= query)
  table(sce$C_SNNGraph.k200, sce$cell_type)
  #      E   I1   I2    M   TC
  # 1    0 1830   31    0  512
  # 2    0   19  762   25  211
  # 3  577    0    0    0  175
  # 4    0    0    8  960  253
  best.matched.clusterID['C_SNNGraph.k200','TC'] = 2  # manually verify
  
  save(best.matched.clusterID, file='best.matched.clusterID.RData')
  best.matched.clusterID
                                                                                           
## 2.3) calculate Jaccard-similarity score for CT identifications  ##
## between each predictions and that of BioTIP using the author-given true_label         ##
## after finding the reference cluster ID (i.e. best-matched IDs) for the QuanTC clusters     ##
## And maually adding the 'best matches' for additional BioTIP/QuanTC runs of this datasets   ##
  jaccard.CT <- array()
  n=2  # assuming TC and M are the transition state !!!!!!
  for(i in 1:length(CT)){
    x <- names(CT)[i]
    reference <- best.matched.clusterID[x,1:n]
    jaccard.CT[i] <- jaccard.sim(reference, CT[[x]]) 
  }
  names(jaccard.CT) <- names(CT)
  round(jaccard.CT, 2)
  # C_cell_type    C_Leiden_0.4    C_Leiden_0.8    C_Leiden_1.2 
  # 0.5             0.5             0.0             0.0 
  # C_SNNGraph.k100  C_SNNGraph.k20 C_SNNGraph.k200      QuanTC_run 
  # 0.5             0.0             0.5             1 

  
 ## 3.4) evaluate the F1 score for the identified CTS genes    ##
 ## for QuanTC identification, we ran k=4 and got only one transition genes (Sha 2020, Fig 2a)
 
  ## we first generate a proxy 'golden standard' of CTS genes,                                
  ## which were 10 genes identifed by at least 4 out of 4 BioTIP's identfications excluding QuanTC -------------------------- 
  {  
    CTS.reference <- read.table('QuanTC_run/Fig2a_genes.txt')[,1]
    CTS.reference
    #[1] "snailt"  "SNAIL"   "miR34t"  "zebt"    "miR200t" "tgfR"    "tgft"    "Ecad"    "Ovol2"
    
    n.select=4
    select <- c( "cell_type", "Leiden_0.4",  "SNNGraph.k100","SNNGraph.k200")  #,"QuanTC_run")
    GS.CTS <- NULL
    for(m in 1:length(select)){
       if(grepl('QuanTC', select[m])) {
          set2 <- CTS.reference 
          #set2 <- read.table('QuanTC_run/transition_gene_I1.I2.txt')[,1]
          GS.CTS = c(GS.CTS, set2)
        } else  {
          load(file=paste0('C_',select[m],'/BioTIP.res.RData'))
          set2 <- res$CTS.candidate[which(res$significant)]
          best.match <- best.matched.clusterID[paste0('C_',select[m]), ]   
          x <- grep(best.match, names(set2))
          GS.CTS = c(GS.CTS, unlist(set2[x]))
        }
      }
      x <- table(GS.CTS)
      GS.CTS <- names(x)[which(x>= n.select)]
      GS.CTS
      #  [1] "miR200t" "tgfR"    "tgft"    "ZEB"     "zebt"    "ZR1"     "ZR2"     "ZR3"     "ZR4"    
      # [10] "ZR5"
 
     (n.GS <- length(GS.CTS))  # 4
      write.table(GS.CTS, file=paste0('GS.CTS_at.least.',n.select,'identificaitons_fr.',length(select),'predicitons.txt'),
                  row.names=F, col.names=F)
  }
  

  # load the CTS prediction with different clustering methods
  running.methods <- list.dirs()
  running.methods <- running.methods[-1]
  methods <- lapply(running.methods, function(x) unlist(strsplit(x, "./C_"))[2]) %>% unlist()
  x <- which(is.na(methods))
  if(length(x)>0) { 
    methods <- methods[-x]
    running.methods <- running.methods[-x]
  }
  ## manually adding the QuanTC's method
  running.methods <- c(running.methods, './QuanTC_run')
  methods <- c(methods,'QuanTC_run')
  
  res.TC <- get.multiple.F1.score(match.col.name=c('TC'),
                               CTS.reference = GS.CTS,
                               QuanTC.genes.list.id = 1, ## using the published 9 gene (Fig for QuanTC results
                               running.methods=running.methods, 
                               methods=methods,
                               sce=sce)
  res.I1 <- get.multiple.F1.score(match.col.name=c('I1'),
                                  CTS.reference = GS.CTS,
                                  QuanTC.genes.list.id = 1,
                                  running.methods=running.methods, 
                                  methods=methods,
                                  sce=sce)
  
  ## 3.5) reproting table #######################
  tb <- data.frame(Jaccard.CT = round(jaccard.CT,3), 
                   F1.TC = round(res.TC$F1.scores$TC, 3), 
                   F1.TC.ctl = round(res.TC$F1.ctl$TC,3),
                   F1.I1 = round(res.I1$F1.scores$I1, 3), 
                   F1.I1.ctl = round(res.I1$F1.ctl$I1,3)
                   )
  nrow(tb) #[1] 8
  tb$detected.CT <- lapply(CT, toString) %>% unlist()
  tb$best.matched.cluster <- apply(best.matched.clusterID,1, toString)

  Normalized.F1 <- Normalize.F1(c(tb$F1.TC, tb$F1.TC.ctl,
                                  tb$F1.I1, tb$F1.I1.ctl))
  tb$Norm.F1.TC = Normalized.F1[1:8]
  tb$Norm.F1.TC.ctl = Normalized.F1[9:16]
  tb$Norm.F1.I1 = Normalized.F1[17:24]
  tb$Norm.F1.I1.ctl = Normalized.F1[25:32]
  
  (mytable <- tb[,
                 c('best.matched.cluster','detected.CT','Jaccard.CT', 
                   'F1.TC.ctl', 'F1.TC','F1.I1.ctl', 'F1.I1',
                   'Norm.F1.TC.ctl', 'Norm.F1.TC','Norm.F1.I1.ctl', 'Norm.F1.I1')])
  #                   best.matched.cluster detected.CT Jaccard.CT F1.TC.ctl F1.TC F1.I1.ctl F1.I1
  # C_cell_type                     2, 1          TC        0.5     0.353 0.769     0.000 0.000
  # C_Leiden_0.4                    3, 1           4        0.5     0.438 0.769     0.000 0.000
  # C_Leiden_0.8                    4, 3                    0.0     0.000 0.000     0.000 0.000
  # C_Leiden_1.2                    4, 1          14        0.0     0.000 0.000     0.000 0.000
  # C_SNNGraph.k100                 2, 3           3        0.5     0.353 0.769     0.000 0.000
  # C_SNNGraph.k20                  7, 2                    0.0     0.000 0.000     0.000 0.000
  # C_SNNGraph.k200               TC, I1           2        0.5     0.471 0.667     0.000 0.000
  # QuanTC_run                    TC, I1          I1        0.5     0.583 0.267     0.267 0.267
  
  write.table(mytable,   file='jaccard.CT_F1.CTS.txt', 
              sep='\t')


}
###############
## PLOT #######
###############
setwd('F:/projects/BioTIP/result/QuanTC_simulation/simulatedEMT_robustness/')
## generate bar plot #######################
df <- read.table('jaccard.CT_F1.CTS.txt', header=T, sep='\t')
df$method = rownames(df)
dim(df)
#[1] 8 12
df$method[1:7] =lapply(rownames(df)[1:7], function(x) unlist(strsplit(x,split='C_',fixed=TRUE))[2]) %>% unlist()
df$method=factor(df$method, levels=df$method)

p1 <- ggbarplot(df, "method", "Jaccard.CT", orientation = "horiz") +
  geom_col(aes(fill = Jaccard.CT)) + 
  scale_fill_gradient2(low = "white",high = "green") +
  theme_bw(base_size=10)  +        # use a black-and-white theme with set font size
  theme(legend.position = "none") + coord_flip()   
p2 <- ggbarplot(df, "method", "Norm.F1.TC", orientation = "horiz") +
  geom_col(aes(fill = Norm.F1.TC)) + 
  scale_fill_gradient2(low = "white",high = "blue") +
  theme_bw(base_size=10)    +      # use a black-and-white theme with set font size
  theme(legend.position = "none")  + coord_flip()  


pdf(file='bar_robustness.pdf', width=6, height=4)
gridExtra::grid.arrange(
  p1, p2#, p3  
  ,ncol=2)
dev.off()




{
df <- data.frame(F1 = c(unlist(res.TC$F1.ctl),unlist(res.TC$F1.scores),
                        unlist(res.I1$F1.ctl),unlist(res.I1$F1.scores) ),
                 CT = c(rep('TC', sum(lengths(res.TC$F1.ctl))+sum(lengths(res.TC$F1.scores))),
                        rep('I1', sum(lengths(res.I1$F1.ctl))+sum(lengths(res.I1$F1.scores))) ),
                 CTS = c(rep('ctl', sum(lengths(res.TC$F1.ctl))), rep('CTS', sum(lengths(res.TC$F1.scores))),
                         rep('ctl', sum(lengths(res.I1$F1.ctl))), rep('CTS', sum(lengths(res.I1$F1.scores))) ) 
)
df$CTS <- factor(df$CTS, levels = c('ctl','CTS')) 
df$CT <- factor(df$CT, levels =c('TC','I1'))

df$Norm.F1 = Normalize.F1(df$F1)
df$color.lab = paste0(df$CT,df$CTS)
# To ensure that easy and difficult CT state have equal influence on the final score, 
# we normalized the scores at each CT state across the different MinModuleSize setting.

library(ggbeeswarm)
library(grid)
library(gridExtra)


  p1 <- ggplot(data=subset(df, F1>0), aes(x=CT, y=F1, color=color.lab)) +
    geom_violin(position = position_dodge(width = 0.5)) +
    geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
    xlab('predicted for cluster') + ylab('F1 score for the CTS') +
    scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue", "grey","red")) 
  
  p2 <- ggplot(data=subset(df, F1>0), aes(x=CT, y=Norm.F1, color=color.lab)) +
    geom_violin(position = position_dodge(width = 0.5)) +
    geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
    xlab('predicted for cluster') + ylab('Normalized F1 score for the CTS') +
    scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue", "grey","red")) 
  
  pdf(file=paste0('CTS.F1.positive_GS_2CTs.pdf'))
  gridExtra::grid.arrange(
    p1, p2, top = "CTS stability", bottom = "different ModuleSize settings"
    ,ncol=1)
  dev.off()


dim(sce)  # 18 5363
## prsent the statistics using proxy gold standard here an in the figure
# ns: p > 0.05
# *: p <= 0.05
# **: p <= 0.01
# ***: p <= 0.001
# ****: p <= 0.0001

tmp = subset(df, CT=='TC' & F1>0)
t.test(tmp$F1[which(tmp$CTS=='CTS')], tmp$F1[which(tmp$CTS=='ctl')])
#t = 0.83669, df = 11.99,  p-value = 0.0387

tmp = subset(df, CT=='I1'& F1>0)
t.test(tmp$F1[which(tmp$CTS=='CTS')], tmp$F1[which(tmp$CTS=='ctl')])
# not enough 'x' observations
}


