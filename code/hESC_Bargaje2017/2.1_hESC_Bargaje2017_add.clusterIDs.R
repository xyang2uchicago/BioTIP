# This code doing:
# 1) Extract cell clusters of the 96 cells using 
# the soft-thresholding clusters with or without the transition cells (by QUANTC), 
## consensus clusters (by SC3), SNNGraph with two parameters (by scran), 
## and Leiden clusters with 4 resolution settings (by Seurat4).
# 2) Save these cluster IDs into the SingleCellExperimentobject.
## and plot umap for each clustering methods
#  
# last update 1/25/2022
# by Holly yang


setwd('F:/projects/BioTIP/result/Bargaje2017_EOMES/')
subDir = 'hESC_Bargaje2017_robustness'

# normalized data were downloaded and reanalyzed
# refer to GSE87038_E8.25_normalize.R

library(dplyr)
library(scater)
library(scran)
library(SC3)
library(Seurat) #Seurat version 4.0.6
#library(leiden) 
#library(monocle3)

parameters = list()
parameters$k = 10 # An integer number of nearest neighboring cells to use when creating the k nearest neighbor graph for Louvain/Leiden/SNNGraph clustering.

################################################################################
## 0) load the R object                                                       ##
## refer to data_matrix_timepoint.R; data_matrix_Cluster_CMpath.R             ##
## focusing on 929 cells analyzed by QuanTC                                   ##
## refer to SC3_GSE87038_E8.25_HEP.R                                          ##
################################################################################

load('F:/projects/BioTIP/result/Bargaje2017_EOMES/sce.RData')
sce
# class: SingleCellExperiment 
# dim: 96 1896 
# metadata(0):
#   assays(1): logcounts
# rownames(96): ACVR1B ACVR2A ... WNT5A WNT5B
# rowData names(0):
#   colnames(1896): Cell1 Cell109 ... Cell2318 Cell2321
# colData names(4): CollectionTime Consensus.Cluster Phenotype.FACS mesodermPath
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):
table(colData(sce)$Consensus.Cluster)
#  1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
# 77 107 130 299 157  64  70  81  67  77 161  86  86  80  37 173  95  49 
table(colData(sce)$CollectionTime)

int <- which(colData(sce)$CollectionTime %in% c('Day0', 'Day1', 'Day1.5', 'Day2', 'Day2.5') |
               (colData(sce)$CollectionTime =='Day3' & colData(sce)$Consensus.Cluster ==10))
sce <- sce[,int]
sce
# class: SingleCellExperiment 
# dim: 96 929 
# metadata(0):
#   assays(1): logcounts
# rownames(96): ACVR1B ACVR2A ... WNT5A WNT5B
# rowData names(0):
#   colnames(929): Cell1 Cell109 ... Cell1421 Cell1429
# colData names(3): CollectionTime Consensus.Cluster Phenotype.FACS
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):
table(colData(sce)$Consensus.Cluster)
#  1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
# 77 107   0   0   0   0   0   0   0   0 161  86  86  80  15 173  95  49
table(colData(sce)$CollectionTime)
#  Day0   Day1 Day1.5   Day2 Day2.5   Day3   Day4   Day5 
#  231    166     93    211    122    106      0      0 


# samplesL <- split(rownames(colData(sce)),f = colData(sce)$Consensus.Cluster)  
# lengths(samplesL)
# #  1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
# # 77 107   0   0   0   0   0   0   0   0 161  86  86  80   15 173  95  49 
# samplesL <- samplesL[-which(lengths(samplesL)< 20)]
# lengths(samplesL)
# #  1  10   2   3   4   5   7   8   9 
# # 77 107 161  86  86  80 173  95  49

#######################################################################################
# 1.1) # extract SNNGraph clusters (by scran) with two settings for the parameter k
# The parameter k indicates the number of nearest neighbors to consider during graph construction, whicc
# we set to 5 for small number of cells (e.g., <500) and 10 to large number of sequenced cells (e.g., 4k).
# In this example, 'PCA' was estimated by Holly with package scater
# by.cluster <- aggregateAcrossCells(sce, ids=colData(sce)$Consensus.Cluster,
#                                   use.assay.type='logcounts', statistics = "sum")
# centroids <- reducedDim(by.cluster, "PCA")
# refer to data_matrix_timpoint.R
########################################################################################

# k: An integer scalar specifying the number of nearest neighbors to consider during graph construction.
SNNGraph.ID <- scran::buildSNNGraph(sce, k= parameters$k, use.dimred = 'PCA')
SNNGraph.ID <- igraph::cluster_walktrap(SNNGraph.ID)$membership

# check the agreeement between new and original clusters using the SNNGraph method
table(as.vector(sce$Consensus.Cluster), SNNGraph.ID)
# SNNGraph.ID
#      1   2   3   4   5   6   7
# 1    1   0   0   0   0   2  74
# 10   0   0   0   0 107   0   0
# 2    2   0   0   0   1   3 155
# 3    0   0   0  71   0  13   2
# 4    3   0   0   2   0  80   1
# 5   73   2   5   0   0   0   0
# 6    2   0   8   5   0   0   0
# 7    0   0 173   0   0   0   0
# 8    0  84  11   0   0   0   0
# 9   40   0   9   0   0   0   0

colData(sce)$C_SNNGraph_k10 = SNNGraph.ID


SNNGraph.ID <- scran::buildSNNGraph(sce, k= 20, use.dimred = 'PCA')
SNNGraph.ID <- igraph::cluster_walktrap(SNNGraph.ID)$membership

# check the agreeement between new and original clusters using the SNNGraph method
table(as.vector(sce$Consensus.Cluster), SNNGraph.ID)
# SNNGraph.ID
#      1   2   3   4   5
# 1    0   2   0   0  75
# 10   0   0   0 107   0
# 2    0   5   0   1 155
# 3    0  15  70   0   1
# 4    0  82   2   0   2
# 5   13  67   0   0   0
# 6    7   5   3   0   0
# 7  173   0   0   0   0
# 8   95   0   0   0   0
# 9   24  25   0   0   0

colData(sce)$C_SNNGraph_k20 = SNNGraph.ID


#save(sce, file=paste0(subDir,'/sce_hESC.RData'), compress=TRUE) # !!!!!!!!!!!!!!!!!

################################################################
# 1.2 ) # extract the QUANTC-assigned 4 clusters for 929 cells
# k=4 was optimalized by QuanTC pipeline 
# refer to QuanTC_soft.thresholding_clusters.m
##################################################################
C_TC <- read.table('F:/projects/BioTIP/result/Bargaje2017_EOMES/QuanTC-modified/Output/k4_4_time/C_TC.txt')
C_TC <- C_TC[,1]
length(C_TC) #[1] 929

index_TC <- read.table('F:/projects/BioTIP/result/Bargaje2017_EOMES/QuanTC-modified/Output/k4_4_time/index_TC.txt')
index_TC <- index_TC[,1]
unique(C_TC[index_TC]) # 5 verified the C_TC is the cluster ID generated by QuanTC

## replace the QuanTC.cluster IDs, TC is the last, to be consistented with those shoing in Fig S2
tmp <- data.frame(C_TC)
tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 1, 'C1'))
tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 2, 'C2'))
tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 3, 'C3'))
tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 4, 'C4'))
tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 5, 'TC'))

table(as.vector(sce$CollectionTime), tmp[,1])
#         C1  C2  C3  C4  TC
# Day0   230   0   0   0   1
# Day1     1   0 163   0   2
# Day1.5   0   1  65   2  25
# Day2     0 114  12  47  38
# Day2.5   0  38   9  55  20
# Day3     0   0   0   0 106


colData(sce)$C_Soft <- tmp[,1]

#save(sce, file=paste0(subDir,'/sce_hESC.RData'), compress=TRUE) # !!!!!!!!!!!!!!!!!


################################################################
# 1.3) # extract consensus clusters  
# refer to http://127.0.0.1:11637/library/SC3/doc/SC3.html
# needs to customize ks accordingily per dataset, the larger range the longer running time.
# In this case, ks=4:10 are tested, 
# and for the related soft-thresholding clustering (QuanTC method), 
# we had take average of the Consensus.Cluster-agreeable clustering results of k=4:10 to get cell-cell similarity matrix M
##################################################################

# # BEGIN NOT repeat, runs hours for big datasets !!!! 
# refer to F:\projects\BioTIP\source\Bargaje2017_EOMES\SC3_cluster.R
# define feature names in feature_symbol column
sce.sc3 = sce
rowData(sce.sc3)$feature_symbol <- rownames(sce.sc3)
# remove features with duplicated names
any(duplicated(rowData(sce.sc3)$feature_symbol))  # F
range(logcounts(sce.sc3))
# [1]     0.00 14.92

### to run SC3 successfully, transform sparsematrix to matrix  !!!!!!!!!!!!!
logcounts(sce.sc3) <- as.matrix(logcounts(sce.sc3))

sum(logcounts(sce.sc3)<1e-16,2)/nrow(logcounts(sce.sc3))>0.95  # TRUE therefore no cell being filtered


# NOT repeat, runs hours !!!! 
# # biology: boolean parameter, defines whether to compute differentially expressed genes, marker genes and cell outliers.
set.seed(2020)
sce.sc3 <- sc3(sce.sc3, ks = 3:10, biology = FALSE) #, svm_max = 5000 is default
# Setting SC3 parameters...
# Your dataset contains more than 2000 cells. Adjusting the nstart parameter of kmeans to 50 for faster performance...
# Calculating distances between the cells...
# Performing transformations and calculating eigenvectors...
# Performing k-means clustering...
# Calculating consensus matrix...

traceback()

# When the sce.sc3 object is prepared for clustering, SC3 can also estimate the optimal number of clusters k in the dataset
# NOT repeat, runs 10 mins  !!!!
sce.sc3 <- sc3_estimate_k(sce.sc3)
str(metadata(sce.sc3)$sc3)
# $ k_estimation   : num 32

# to save space, transform back matrix to sparse matrix
counts(sce.sc3) <- as(counts(sce.sc3), 'dgCMatrix')
logcounts(sce.sc3) <- as(logcounts(sce.sc3), 'dgCMatrix')
# 
# sce <- sce.sc3
# save(sce, file='F:/projects/BioTIP/result/Bargaje2017_EOMES/Cluster/sce_SC3.RData') ##!!!!!!!!!!!!!!!!!!!
# gc()
# END DO NOT REPET !!!!!!!!!!!!!

#load('F:/projects/scRNA/results/GSE87038_gastrulation/QuanTC/QuanTC_HEP/sce_E8.25HEP_uncorrected_SC3.RData')
# sce.sc3 <- sce
# load(file='hESC_Bargaje2017_robustness/sce_E8.25_HEP.RData')
# colData(sce) = colData(sce)[,1:22]

## manually pick the optimal matches to follow up
table(as.vector(sce$Consensus.Cluster), colData(sce.sc3)$sc3_6_clusters)
#      1   2   3   4   5   6
# 1   75   0   0   2   0   0
# 10   0   0   0   0   0 107
# 2  154   0   3   3   0   1
# 3    3   0  77   6   0   0
# 4    0   0  84   2   0   0
# 5    0   6   0  69   5   0
# 6    0  10   4   1   0   0
# 7    0 157   0   0  16   0
# 8    0   0   0   0  95   0
# 9    0  20   0  28   1   0


# load(file='sce_E8.25_HEP.RData')
colData(sce)$C_consensus_ks4 = colData(sce.sc3)$sc3_4_clusters
colData(sce)$C_consensus_ks5 = colData(sce.sc3)$sc3_5_clusters
colData(sce)$C_consensus_ks6 = colData(sce.sc3)$sc3_6_clusters
colData(sce)$C_consensus_ks7 = colData(sce.sc3)$sc3_7_clusters
colData(sce)$C_consensus_ks8 = colData(sce.sc3)$sc3_8_clusters
colData(sce)$C_consensus_ks9 = colData(sce.sc3)$sc3_9_clusters
colData(sce)$C_consensus_ks10 = colData(sce.sc3)$sc3_10_clusters


rm(sce.sc3)

# save(sce, file=paste0(subDir,'/sce_hESC.RData'), compress=TRUE) # !!!!!!!!!!!!!!!!!


################################################################
# 1.4.1) # extract Leiden clustering (using Seurat)  
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
# generate a pseudo count to run Seurat
counts(sce) <- as.matrix(2^logcounts(sce)-1)

# convert from SingleCellExperiment
sce.seurat <- as.Seurat(sce)
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
sce.seurat
# An object of class Seurat 
# 96 features across 929 samples within 1 assay 
# Active assay: originalexp (96 features, 0 variable features)
# 3 dimensional reductions calculated: PCA, TSNE, UMAP


# Computes the k.param nearest neighbors for a given dataset
sce.seurat <- FindNeighbors(sce.seurat, reduction = "PCA", k.param = parameters$k)

sce.seurat <- FindClusters(sce.seurat, resolution = 0.4, algorithm = 4) # smaller  number of communities 
table(Idents(sce.seurat), as.vector(colData(sce)$Consensus.Cluster))  
#    1  10   2   3   4   5   6   7   8   9
# 1  73   0 151   1   0   0   0   0   0   0
# 2   0   0   0   0   0   2  10 172   1   2
# 3   1   0   2   0   1  75   0   1   0  47
# 4   0 107   1   0   0   0   0   0   0   0
# 5   3   0   7  13  81   0   2   0   0   0
# 6   0   0   0   0   0   3   0   0  94   0
# 7   0   0   0  72   4   0   3   0   0   0
colData(sce)$C_Leiden_0.4 = Idents(sce.seurat)

sce.seurat <- FindClusters(sce.seurat, resolution = 0.8, algorithm = 4) # smaller  number of communities 
table(Idents(sce.seurat), as.vector(colData(sce)$Consensus.Cluster))  
#     1  10   2   3   4   5   6   7   8   9
# 1    0   0   0   0   0   2   0 164   0   2
# 2   18   0 122   1   0   0   0   0   0   0
# 3    0 107   1   0   0   0   0   0   0   0
# 4    1   0   5  13  81   0   2   0   0   0
# 5    0   0   0   0   0   3   0   0  94   0
# 6    1   0   2   0   1  67   0   0   0  24
# 7   57   0  31   0   0   0   0   0   0   0
# 8    0   0   0  72   4   0   3   0   0   0
# 9    0   0   0   0   0   8   0   1   0  23
# 10   0   0   0   0   0   0  10   8   1   0

colData(sce)$C_Leiden_0.8 = Idents(sce.seurat)

# sce.seurat <- FindClusters(sce.seurat, resolution = 0.6, algorithm = 4) # smaller  number of communities 
# table(Idents(sce.seurat), as.vector(colData(sce)$C_Leiden_0.8))  # the same
sce.seurat <- FindClusters(sce.seurat, resolution = 1.2, algorithm = 4) # smaller  number of communities 
table(Idents(sce.seurat), as.vector(colData(sce)$C_Leiden_0.8))  # the same

table(Idents(sce.seurat), as.vector(colData(sce)$Consensus.Cluster))  
#     1  10   2   3   4   5   6   7   8   9
# 1    0   0   0   0   0   2   0 142   0   2
# 2    0 107   1   0   0   0   0   0   0   0
# 3    1   0   5  13  81   0   2   0   0   0
# 4    1   0   2   0   1  67   0   0   0  24
# 5   11   0  77   0   0   0   0   0   0   0
# 6   56   0  32   0   0   0   0   0   0   0
# 7    0   0   0  72   4   0   3   0   0   0
# 8    0   0   0   0   0   3   0   0  66   0
# 9    8   0  44   1   0   0   0   0   0   0
# 10   0   0   0   0   0   8   0   1   0  23
# 11   0   0   0   0   0   0   0   0  28   0
# 12   0   0   0   0   0   0   0  22   0   0
# 13   0   0   0   0   0   0  10   8   1   0
colData(sce)$C_Leiden_1.2 = Idents(sce.seurat)

# In this case, remove teh pseudo counts
counts(sce) <- NULL
#save(sce, file=paste0(subDir,'/sce_hESC.RData'), compress=TRUE) # !!!!!!!!!!!!!!!!!



####################################################################
## 2) Save these cluster IDs into the SingleCellExperimentobject. ##
## plot different clustering restuls                              ##
####################################################################
save(sce, file=paste0(subDir,'/sce_hESC.RData'), compress=TRUE) # !!!!!!!!!!!!!!!!!

rm(cds)


x <- grep('C_', colnames(colData(sce)))
(n=length(x)) # 13

pdf(file=paste0(subDir,"/TSNE.",subDir,"_clustering_methods.pdf"), width=10, height=9)
gridExtra::grid.arrange(
   plotReducedDim(sce, dimred='TSNE', colour_by='CollectionTime', #add_legend=FALSE,
                 text_by='CollectionTime', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('CollectionTime')  #+ ylim(5,20)  
   ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Soft', #add_legend=FALSE,
                  text_by='C_Soft', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('Soft clustering by QuanTC k4') #+ ylim(5,20)
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_SNNGraph_k10', #add_legend=FALSE,
                  text_by='C_SNNGraph_k10', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('SNNGraph_k10')  #+ ylim(5,20)
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_SNNGraph_k20', #add_legend=FALSE,
                  text_by='C_SNNGraph_k10', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('SNNGraph_k20')  #+ ylim(5,20)
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks4', #add_legend=FALSE,
                  text_by='C_consensus_ks4', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('consensus_ks4') #+ ylim(5,20) 
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks9', #add_legend=FALSE,
                  text_by='C_consensus_ks9', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('consensus_ks9') #+ ylim(5,20)
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Leiden_0.4', #add_legend=FALSE,
                  text_by='C_Leiden_0.4', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('Leiden_0.4') #+ ylim(5,20)
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Leiden_0.8', #add_legend=FALSE,
                  text_by='C_Leiden_0.8', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('Leiden_0.8') #+ ylim(5,20)
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Leiden_1.2', #add_legend=FALSE,
                  text_by='C_Leiden_1.2', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('Leiden_1.2') #+ ylim(5,20)
  ,ncol=3
)
gridExtra::grid.arrange(
  plotReducedDim(sce, dimred='TSNE', colour_by='Consensus.Cluster', #add_legend=FALSE,
               text_by='Consensus.Cluster', text_size = 4, text_colour='black', point_size=0.5) +
  ggtitle('CollectionTime') 
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks5', #add_legend=FALSE,
                  text_by='C_consensus_ks5', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('consensus_ks5') #+ ylim(5,20) 
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks6', #add_legend=FALSE,
                  text_by='C_consensus_ks6', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('consensus_ks6') #+ ylim(5,20) 
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks7', #add_legend=FALSE,
                  text_by='C_consensus_ks7', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('consensus_ks7') #+ ylim(5,20) 
  ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks8', #add_legend=FALSE,
                  text_by='C_consensus_ks8', text_size = 4, text_colour='black', point_size=0.5) + 
    ggtitle('consensus_ks8') #+ ylim(5,20) 
  ,ncol=3, nrow=3)

dev.off()




