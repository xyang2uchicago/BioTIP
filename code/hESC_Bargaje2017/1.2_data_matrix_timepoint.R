# R4.0.2
 setwd('F:/projects/BioTIP/result/Bargaje2017_EOMES/timepoint/')

 library(RColorBrewer)
 library(SingleCellExperiment)
 library(org.Mm.eg.db)
 library(scater)
 library(scran)
 library(pheatmap)
 library(limma)
 library(edgeR)
 library(dynamicTreeCut)
 library(cluster)
 library(bluster)
 
 library(MouseGastrulationData)
 head(EmbryoCelltypeColours)
 length(EmbryoCelltypeColours) # 37
 
###################################################
## 1) read the published data matrix of 96 genes ##
###################################################
# Dataset S2: This file contains all necessary information to re-analyze the data. 
# We provide the metadata and log2Ex values for all cells 
# (after removing 38 cells for quality control) 
# (processed as described in Materials and Methods). 
# Phenotypes.FACS = phenotypic state based on flow cytometry measurements (Fig. S1)																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																											
 dat <-  read.table('F:/projects/BioTIP/ref/development/Bargaje2017PNAS_SData_logExp.txt',sep="\t",header=T)
 cli <- t(data.frame(dat[1:3,]))
 cli <- cli[-1,]
 dim(cli)  #[1] 1896    3
 head(cli)
 #      CollectionTime Consensus Cluster Phenotype.FACS
 # Cell1   "Day0"         "1"               "d0.cKIT-h"   
 # Cell109 "Day0"         "1"               "d0.meso-"    
 # Cell11  "Day0"         "1"               "d0.cKIT-h"   
 # Cell111 "Day0"         "1"               "d0.meso-"    
 # Cell112 "Day0"         "1"               "d0.meso-"    
 # Cell114 "Day0"         "1"               "d0.meso-"
 colnames(cli)[2] <- 'Consensus.Cluster'
 cli <- as.data.frame(cli)
 cli$CollectionTime <- as.factor(cli$CollectionTime)
 
 table(cli[,1])
 # Day0   Day1 Day1.5   Day2 Day2.5   Day3   Day4   Day5 
 # 231    166     93    211    122    258    456    359 
 table(cli[,2])
 #  1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
 # 77 107 130 299 157  64  70  81  67  77 161  86  86  80  37 173  95  49 
table(cli[,3])
 # d0.cd13+     d0.cKIT-h     d0.cKIT-l     d0.cKIT-m      d0.meso- 
 #   35            30            20            22            39 
 # d0.Tra1-60-   d0.Tra1-60+    d1.5.cd13+   d1.5.cKIT-h   d1.5.cKIT-l 
 # 41            44            17            25             9 
 # d1.5.cKIT-m    d1.5.meso-    d1.5.meso+      d1.cd13+     d1.cKIT-h 
 # 30             7             5            27            21 
 # d1.cKIT-l     d1.cKIT-m      d1.meso-   d1.Tra1-60-   d1.Tra1-60+ 
 #   11            23            30            24            30 
 # d2.5.cd13+   d2.5.cKIT-h   d2.5.cKIT-l   d2.5.cKIT-m    d2.5.meso- 
 #   11            32            17            41            10 
 # d2.5.meso+     d2.cKIT-h     d2.cKIT-l     d2.cKIT-m       d2.KDR- 
 #   11            25            21            26            36 
 # d2.KDR+      d2.meso+      d2.ror2+  d3.cd13.kdr+ d3.cd13h.kdr- 
 #   37            38            28            92            77 
 # d3.cd13l.kdr-     d4.Hybrid       d4.kdr- d4.kdr-cxcr4+ d4.kdr+cxcr4- 
 #   89           244            74            64            74 
 # d5.Hybrid d5.kdr.cxcr4- d5.kdr.cxcr4+ d5.kdr+cxcr4- 
 #   146            76            73            64 

dat <- dat[-1:4,]
rownames(dat) <- dat[,1]
dat <- dat[,-1]
sce <- SingleCellExperiment(assays=SimpleList(logcounts=dat),
                            colData=cli)
colData(sce) <- cli
dim(sce) # 96 1896
rm(cli, dat)



###########################################
## 2) reduce dimention and visualization ##
###########################################
set.seed(100) # See below.
sce <- runPCA(sce) 
# You're computing too large a percentage of total singular values, use a standard svd instead.
sce <- runTSNE(sce, dimred="PCA")
#plotReducedDim(sce, dimred="TSNE", colour_by=clust)
sce <- runUMAP(sce, dimred="PCA" )
#plotReducedDim(sce, dimred="UMAP", colour_by=clust)


pdf(file='ReduceDim.pdf', height=4.5, width=10) #!!!!!!!!!!!!!!!!!

MyEmbryoCelltypeColours <- EmbryoCelltypeColours[1:length(unique(sce$Consensus.Cluster))]
names(MyEmbryoCelltypeColours) <- 1:length(unique(sce$Consensus.Cluster))

# Performing the calculations on the PC coordinates, like before.
sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=sce$Consensus.Cluster)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, sce$Consensus.Cluster, sil.data$other))
sil.data$cluster <- factor(sce$Consensus.Cluster)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="PCA", colour_by='CollectionTime', text_by='Consensus.Cluster',point_size=0.5),# +
    #scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "TSNE"), clusters=sce$Consensus.Cluster)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, sce$Consensus.Cluster, sil.data$other))
sil.data$cluster <- factor(sce$Consensus.Cluster)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="TSNE", colour_by='CollectionTime', text_by='Consensus.Cluster',point_size=0.5), #+
    #scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "UMAP"), clusters=sce$Consensus.Cluster)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, sce$Consensus.Cluster, sil.data$other))
sil.data$cluster <- factor(sce$Consensus.Cluster)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="UMAP", colour_by='CollectionTime', text_by='Consensus.Cluster',point_size=0.5),# +
   # scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ncol=2
)

dev.off()

colData(sce)$Consensus.Cluster <- factor(colData(sce)$Consensus.Cluster,
                                         levels = c(1:10,11:18))

################################################################################
## 3) defining cell identity by the expression levels of known lineage-marker ##
################################################################################
newmarkers = c("NANOG","TBX2", "GATA4","MESP1","GSC","KDR","HAND1","SOX17","EOMES","T")

plotExpression(sce, 
               features=newmarkers, 
               x="Consensus.Cluster", colour_by="Consensus.Cluster",
               point_size=0.05)
dev.copy2pdf(file='../Marker_expression.pdf')


###################################################
## 4) trajectory analysis, Disagree with biology ##
###################################################

## Minimum spanning tree constructed trajectory

library(scater)

## by clusters , statistics == mean or mdian didn't change patter !------------------
by.cluster <- aggregateAcrossCells(sce, ids=colData(sce)$Consensus.Cluster,
                                   use.assay.type='logcounts', statistics = "sum")
 
centroids <- reducedDim(by.cluster, "PCA")
dmat <- dist(centroids)
dmat <- as.matrix(dmat)
set.seed(1000)
g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
mst <- igraph::minimum.spanning.tree(g)


pdf(file="trajectory_scater.pdf")
plot(mst)

pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
coords <- reducedDim(by.cluster, "PCA")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)
plotReducedDim(sce, 'PCA', colour_by="Consensus.Cluster", 
               text_by="Consensus.Cluster", text_size = 1) + 
geom_line(data=stuff, mapping=aes(x=PC1, y=PC2, group=group))

pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
coords <- reducedDim(by.cluster, "TSNE")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)
plotReducedDim(sce, 'TSNE', colour_by="Consensus.Cluster", 
               text_by="Consensus.Cluster", text_size = 1) + 
geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))

pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
coords <- reducedDim(by.cluster, "UMAP")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)
plotReducedDim(sce, 'UMAP', colour_by="Consensus.Cluster", 
               text_by="Consensus.Cluster", text_size = 1) + 
  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))

dev.off()

plotReducedDim(sce, dimred='TSNE', colour_by='BMP2', swap_rownames=NULL, ncomponents=c(1,2),
               text_by='Consensus.Cluster', text_size = 2, text_colour='red', point_size=1, 
               by_exprs_values='logcounts') 
 


## begin [if rerun] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
load('../sce.RData')

## end of [if rerun]
myplot<- function(title ='', dimred.sce, dimred = 'TSNE', 
                  newmarkers ='Kdr', 
                  ncomponents = c(1,2), 
                  text_by ='label',
                  swap_rownames ='Symbol', 
                  by_exprs_values = "logcounts", 
                  n = 10, 
                  text_size = 5)
{
   if(length(newmarkers) < n*3) newmarkers <- c(newmarkers, rep(newmarkers[1], n*3-length(newmarkers)))
   b <- ceiling(length(newmarkers)/n)-1
   for(x in 0:b)
   {
    gridExtra::grid.arrange(
      plotReducedDim(dimred.sce, dimred=dimred, colour_by=newmarkers[1+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                     text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[1+n*x])) # + geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred,colour_by=newmarkers[2+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                      text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[2+n*x]))#  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred,colour_by=newmarkers[3+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                      text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[3+n*x]))#  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred,colour_by=newmarkers[4+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                      text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[4+n*x]))#  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred, colour_by=newmarkers[5+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                     text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[5+n*x])) # + geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred,colour_by=newmarkers[6+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                      text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[6+n*x]))#  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred,colour_by=newmarkers[7+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                      text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[7+n*x]))#  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred,colour_by=newmarkers[8+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                      text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[8+n*x]))#  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred,colour_by=newmarkers[9+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                      text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
        ggtitle(paste(title,newmarkers[9+n*x]))#  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
      ,plotReducedDim(dimred.sce, dimred=dimred,colour_by=newmarkers[10+n*x], swap_rownames=swap_rownames,ncomponents=ncomponents,
                      text_by=text_by, text_size = text_size, text_colour='red',point_size=1, by_exprs_values=by_exprs_values) + 
                      ggtitle(paste(title,newmarkers[10+n*x]))  
      ,ncol=5
    )
   }
}
newmarkers = unique(c("NANOG","EOMES","T","GSC","MESP1","GATA4","KDR","HAND1","TBX2", "SOX17",
                      "BMP4","FGF8",  # , "GBX2", "OTC2", "TP73L", "PAX2", "PAX6", "SOX1") #Early Ectodermal Lineage Markers
                      "DKK1",  "TNNT2","WNT3A", "KIT", "TBX5","PTX2","ISL1", 
                      'BMP2',"WNT5A",'FGF8','HAND2','HRT2', 'MESP2','DLL1','DLL3','SFRP1')) # 'HRT2'=='HEY2' 
length(newmarkers)
#[1] 27
setdiff(newmarkers, rownames(sce))
#character(0)
pdf(file="../tSNE.iPSC_marker.pdf", width=14)
myplot(title='', dimred.sce=sce, newmarkers=newmarkers, 
       text_by='Consensus.Cluster', swap_rownames=NULL, by_exprs_values='logcounts')
dev.off()


## MEgan's differentiated markers 
newmarkers = unique(c("BMP4","EOMES",  "WNT3A","TBX20","TBX5",
                      "HAND1","TNNT2","MYL4","NKX2.5","HRT2",
                      "WNT5A", "HAND2","MESP2","SFRP1","DLL1","DLL3",
                      "BMP2", "FGF8", "GATA4", "TBX2", "KDR","SHH",
                      ## following are CTS_C10_19g
                      "BMP2",  "DLL1" , "DLL3",  "EVX1",  "FGF12", "FGF8",  "FOXC1",
                      "FZD7",  "GATA4", "HAND2", "HRT2",  "ISL1",  "MESP1", "MESP2",
                      "MIXL1", "SFRP1", "TBX2",  "WNT4",  "WNT5A"
)) # 'HRT2'=='HEY2' 
pdf(file="../tSNE.iPSC_marker_differentated.pdf", width=14)
myplot(title='', dimred.sce=sce, newmarkers=newmarkers, 
       text_by='Consensus.Cluster', swap_rownames=NULL, by_exprs_values='logcounts')
dev.off()



newmarkers = unique(c('DLL3', 'FZD1', 'FZD2', 'SFRP1',
                        'ALCAM', 'BAMBI', 'EPCAM', 'HRT2', 
                        'NANOG', 'PDGFA', 'PDGFB', 'PDGFRB')
) # 'HRT2'=='HEY2' 
plotExpression(sce, 
               features=newmarkers, 
               x="CollectionTime", colour_by="CollectionTime",
               point_size=0.05)
dev.copy2pdf(file='../Marker_expression_timpoint.pdf')

newmarkers = unique(c('DLL3', 'FZD1', 'FZD2', 'SFRP1',
                      'ALCAM', 'BAMBI', 'EPCAM', 'HRT2', 
                      'NANOG', 'PDGFA', 'PDGFB', 'PDGFRB')
) # 'HRT2'=='HEY2' 
plotExpression(sce[,which(sce$mesodermPath)], 
               features=newmarkers, 
               x="CollectionTime", colour_by="CollectionTime",
               point_size=0.05)
dev.copy2pdf(file='../Marker_expression_timpoint_929cells.pdf')




# library(cowplot)
# library(dplyr)
# library(scater)
# library(Seurat)
# library(cowplot)
# 
# sce.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")
# # Error: No data in provided assay - counts
# 
# newmarkers = c("EOMES","BMP4","WNT3A",'HAND1')
# # Plot the expression of each of the genes of interest on the tSNE
# FeaturePlot(sce, features = newmarkers)


############################
## 4) BioTIP application  ##
############################
library(BioTIP)
packageVersion('BioTIP') # '1.5.0'

# Ic, we generated time point-specific Log2Ex matrices taht were downloaded from teh publication 
# For each of them (days 0, 1, 1.5, 2, 2.5, and 3/M-only mesoderm-specific cells)
load('../sce.RData')
tmp <- subset(colData(sce), CollectionTime=='Day3' )
table(tmp$Consensus.Cluster)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
# 0   0   0   0   0  22   0   0   0 106 130   0   0   0   0   0   0   0 

int <- which(colData(sce)$CollectionTime %in% c('Day0', 'Day1', 'Day1.5', 'Day2', 'Day2.5') |
               (colData(sce)$CollectionTime =='Day3' & colData(sce)$Consensus.Cluster ==10))
length(int)  #[1] 929

table(colData(sce)[int,]$CollectionTime)
#  Day0   Day1 Day1.5   Day2 Day2.5   Day3   Day4   Day5 
#  231    166     93    211    122    106      0      0 

colData(sce)$mesodermPath = FALSE
colData(sce)$mesodermPath[int] = TRUE
all(which(colData(sce)$mesodermPath) == int)  #TRUE
save(sce, file='../sce.RData')  #!!!!!!!!!!!!!!!!!!!!!!!



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
  

# 1) module construction

# Pre-selection Transcript
cut.preselect = 0.8  # out of 96 genes
cut.fdr = 0.2   
cut.minsize = 10  

samplesL <- split(rownames(colData(sce)),f = colData(sce)$CollectionTime)
lengths(samplesL)
# Day0   Day1 Day1.5   Day2 Day2.5   Day3   Day4   Day5 
# 231    166     93    211    122    106      0      0 
samplesL <- samplesL[-which(lengths(samplesL)<10)]
lengths(samplesL)
# Day0   Day1 Day1.5   Day2 Day2.5   Day3 
# 231    166     93    211    122    106 

logmat <- as.matrix(logcounts(sce))
dim(logmat) # [1] 96 929

#this optimizes sd selection
set.seed(2020)
testres <- optimize.sd_selection(logmat, samplesL, B=100, cutoff=cut.preselect,
                                 times=.75, percent=0.8)
class(testres) #[1] "list"
class(testres[[1]])  #[1] "matrix" "array" 
lapply(testres, dim)
# $Day0
# [1]  76 231
# 
# $Day1
# [1]  76 166
# 
# $Day1.5
# [1] 72 93
# 
# $Day2
# [1]  72 211
# 
# $Day2.5
# [1]  75 122
# 
# $Day3
# [1]  75 106


# Network Partition
igraphL <- getNetwork(testres, fdr = cut.fdr)
# Day0:76 nodes
# Day1:76 nodes
# Day1.5:72 nodes
# Day2:72 nodes
# Day2.5:75 nodes
# Day3:73 nodes

cluster <- getCluster_methods(igraphL)

##### plot network #################
####################################
names(igraphL)
#[1]"Day0"   "Day1"   "Day1.5" "Day2"   "Day2.5" "Day3" 

library('igraph')
tmp = igraphL[["Day3"]]
E(tmp)$width <- E(tmp)$weight*3
V(tmp)$community= cluster[["Day3"]]$membership
mark.groups = table(cluster[["Day3"]]$membership)
#mark.groups = mark.groups[which(mark.groups>=10)]
colrs = rainbow(length(mark.groups), alpha = 0.3)
V(tmp)$label <- NA
plot(tmp, vertex.color=colrs[V(tmp)$community], vertex.size = 5,
     mark.groups=cluster[["Day3"]])
legend(1,1, paste0(names(mark.groups),sep=":",mark.groups), text.col=rainbow(length(mark.groups)))

dev.copy2pdf(file='Network.view_Day3.pdf')


# 2.1) Identifying CTS, new MCI score #########
membersL <- getMCI(cluster,testres, adjust.size = FALSE, fun='BioTIP')
names(membersL)
#[1] "members" "MCI"     "sd"      "PCC"     "PCCo"  
save(membersL, file="membersL.RData")

pdf(file=paste0("MCIBar_", cut.preselect,
                "_fdr",cut.fdr,"_minsize",cut.minsize,".pdf"),
    width=10, height=5)
plotBar_MCI(membersL, ylim=c(0,10), minsize = 30)
plotBar_MCI(membersL, ylim=c(0,10), minsize = 20)
plotBar_MCI(membersL, ylim=c(0,10), minsize = cut.minsize)
plotBar_MCI(membersL, ylim=c(0,10), minsize = 5)

dev.off()
#the numbers in the parenthesis: they are total number of clusters, no control of the cell (sample) size per  cluster


# get the statistics using the MCI system
n.state.candidate = 4
topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], 
                   min=cut.minsize, n=n.state.candidate)
names(topMCI)
#[1] "Day2.5" "Day2"   "Day1.5" "Day1" 

# get the state ID and MCI statistics for the two leading MCI scores per state
maxMCIms <- getMaxMCImember(membersL[["members"]],
                            membersL[["MCI"]],min =cut.minsize, n=3)
names(maxMCIms)
#[1] "idx"             "members"         "2topest.members" "3topest.members"

maxMCI = getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
unlist(maxMCI)
#    Day0     Day1   Day1.5     Day2   Day2.5     Day3 
# 2.876082 3.668823 3.693530 3.997960 4.229456 3.666241 

names(maxMCIms[["members"]][names(topMCI)])
#[1] [1] "Day2.5" "Day2"   "Day1.5" "Day1" 

CTS = getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
# Length: 18
# Length: 23
# Length: 37
# Length: 43

## tmp calculates the number of bars within each named state
(tmp = unlist(lapply(maxMCIms[['idx']][names(topMCI)], length)))
# Day2.5   Day2 Day1.5   Day1 
# 3      2      3      2

## here returns all the groups with exactly 2 bars
(whoistop2nd = names(tmp[tmp==2]))
#[1] "Day2" "Day1"

## here returns all the groups with exactly 3 bars
(whoistop3rd = names(tmp[tmp==3]))
#[1] "Day2.5" "Day1.5"

## add the gene members of the 2nd toppest biomodue in the states with exactly 2 bars
if(length(whoistop2nd)>0)  CTS = append(CTS, maxMCIms[["2topest.members"]][whoistop2nd])
names(CTS)
#[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2"   "Day1" 

## add the gene members of the 2nd toppest biomodue in the states with exactly 3 bars
if(length(whoistop3rd)>0)  CTS = append(CTS, maxMCIms[["2topest.members"]][whoistop3rd])  
names(CTS)
#[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2"   "Day1"   "Day2.5" "Day1.5"

## add the gene members of the 3rd toppest biomodue in the states with exactly 3 bars
if(length(whoistop3rd)>0)  CTS = append(CTS, maxMCIms[["3topest.members"]][whoistop3rd])  
names(CTS)
#[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2"   "Day1"   "Day2.5" "Day1.5" "Day2.5" "Day1.5"


CTS[1:4]
# $Day2.5
# [1] "ALCAM"   "BAMBI"   "EMILIN2" "EPCAM"   "EVX1"    "FGF10"   "HRT2"    "ISL1"    "MSX2"   
# [10] "NANOG"   "PDGFA"   "PDGFB"   "PDGFRB"  "T"       "TBX5"    "TNNT2"   "WNT5B"   "WNT4"   
# 
# $Day2
# [1] "ACVRL1" "BMP2"   "DKK1"   "DLL1"   "DLL3"   "FOXC1"  "FZD2"   "GATA4"  "HAND1"  "LEFTY1"
# [11] "MESP1"  "MESP2"  "MSX1"   "MYL3"   "MYOCD"  "NANOG"  "PDGFB"  "PDGFRA" "SFRP1"  "SOX17" 
# [21] "T"      "TBX20"  "WNT5A" 
# 
# $Day1.5
# [1] "ACVR2A"  "ALCAM"   "BAMBI"   "BMP4"    "BMPR1A"  "BMPR2"   "EMILIN2" "ENG"     "FGF12"  
# [10] "FGF8"    "FGFR1"   "FGFR2"   "FOXH1"   "FSTL1"   "FZD4"    "FZD6"    "HEY1"    "HRT2"   
# [19] "LTBP1"   "MYOCD"   "NANOG"   "NOTCH2"  "NOTCH3"  "NUMB"    "PARD3"   "PDGFB"   "PDGFRB" 
# [28] "PTCH1"   "RCOR2"   "RPL35A"  "SIRPA"   "TGFB2"   "TGFBR1"  "TGFBR2"  "TNNT2"   "TUBB"   
# [37] "VEGFA"  
# 
# $Day1
# [1] "ACVR1B" "ACVR2A" "ALCAM"  "ANF"    "BAMBI"  "BMP4"   "BMPR1A" "BMPR2"  "DLL3"   "ENG"   
# [11] "EOMES"  "EPCAM"  "FSTL1"  "FZD1"   "FZD2"   "FZD4"   "FZD6"   "FZD7"   "HEY1"   "HHIP"  
# [21] "HRT2"   "IRX4"   "KDR"    "KIT"    "LTBP1"  "NANOG"  "NOTCH2" "NOTCH3" "NUMB"   "PARD3" 
# [31] "PDGFA"  "PDGFB"  "PDGFRB" "PTX2"   "RCOR2"  "RPL35A" "SERCA"  "SFRP1"  "SIRPA"  "TGFB1" 
# [41] "TGFBR1" "TUBB"   "VEGFA" 

#### extract CTS scores for each biomodule candidate in the following steps: ####
## first to record the max MCI for the n.state.candidate 
maxMCI <- maxMCI[names(CTS)[1:n.state.candidate]]
maxMCI
# Day2.5     Day2   Day1.5     Day1 
# 4.229456 3.997960 3.693530 3.668823 

## then applendix the 2nd highest MCI score (if existing) for the states with exactly 2 bars
if(length(whoistop2nd)>0) maxMCI <- c(maxMCI, getNextMaxStats(membersL[['MCI']], idL=maxMCIms[['idx']], whoistop2nd))
names(maxMCI)
# [1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2"   "Day1"  

## applendix the 2nd highest MCI score (if existing) for the states with exactly 3 bars
if(length(whoistop3rd)>0) maxMCI <- c(maxMCI, getNextMaxStats(membersL[['MCI']], idL=maxMCIms[['idx']], whoistop3rd))
names(maxMCI)
# [1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2"   "Day1"   "Day2.5" "Day1.5" 

## then applendix the 3rd highest MCI score (if existing) for the states with exactly 3 bars
if(length(whoistop3rd)>0) maxMCI <- c(maxMCI, getNextMaxStats(membersL[['MCI']], idL=maxMCIms[['idx']], whoistop3rd, which.next=3))
maxMCI
#   Day2.5      Day2    Day1.5      Day1      Day2      Day1    Day2.5    Day1.5    Day2.5    Day1.5 
# 4.2294562 3.9979599 3.6935297 3.6688229 1.9068669 0.7693816 3.5453475 3.0876365 3.4842829 1.0100347 


## to ensure the same order between maxMCI  and CTS
all(names(CTS) == names(maxMCI)) # TRUE


save(CTS,  maxMCI, whoistop2nd, whoistop3rd, file=paste0("CTS_maxMCI.RData")) # !!!

## estimate significance NOT repeat (takes a while)
C = 1000

simuMCI <- list()
for(i in 1:length(CTS) ){
  n <- length(CTS[[i]])
  simuMCI[[i]] <- simulationMCI(n, samplesL, logmat,  B=C, fun="BioTIP")
}
names(simuMCI) <- names(CTS)
save(simuMCI, P, file=paste0("GenePermutation_",C,"CTS_timepoint", cut.preselect, "_fdr",cut.fdr,"_minsize",cut.minsize,".RData"))


## the sig CTS candidates were c(1:4,7,9)
load(file=paste0("CTS_maxMCI.RData"))
load(file=paste0("GenePermutation_",C,"CTS_timepoint", cut.preselect, "_fdr",
                 cut.fdr,"_minsize",cut.minsize,".RData"))
P <- maxMCI * NA
par(mfrow=c(2,5))
for(i in 1:length(CTS) ){
  n <- length(CTS[[i]])
  P[i] <- length(which(maxMCI[i] <= simuMCI[[i]][names(CTS)[i],]))/C
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las=2, ylim=c(0,4),
                    main=paste(names(CTS)[i], n, "genes",
                               "\n","vs. ",C, "times of gene-permutation"), which2point=names(CTS)[i])
}
dev.copy2pdf(file=paste0("MCI_GenePermutation_1000_CTS_timepoint", cut.preselect, "_fdr",cut.fdr,"_minsize",cut.minsize,".pdf"))

P
 # [1] 0.000 0.000 0.000 0.000 0.268 1.000 0.000 0.116 0.000 1.000
(sig <- which(P<0.01) ) 
#[1] 1 2 3 4 7 9



####### Verifying Tipping Point using IC* score  #################
######  BioTIP score, permutating genes ##############
C= 1000

set.seed(2021)
IC_new <- simuBioTIP_g <- list()  
for(i in 1:length(sig)){
  n <- length(CTS[sig][[i]])
  IC_new[[i]] <- getIc(logmat, samplesL, CTS[sig][[i]], fun="BioTIP", shrink=TRUE)

  simuBioTIP_g[[i]]  <- simulation_Ic(n, samplesL, logmat, B=C,
                                 fun="BioTIP", shrink=TRUE)
}
names(simuBioTIP_g) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
names(IC_new) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')

save(IC_new, simuBioTIP_g,   file="simuBioTIP_g.RData", compress=TRUE)

names(CTS)[sig]
#[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2.5" "Day2.5"  

load(file="simuBioTIP_g.RData")
pdf(file="BioTIP_vsSimulation_6CTS.pdf", height=10)
par(mfrow=c(3,4))
for(i in 1:length(IC_new)){
  plot_Ic_Simulation(IC_new[[i]], simuBioTIP_g[[i]], las = 2, ylab="Ic*",
                     main= names(IC_new)[i],
                     fun="boxplot", which2point= names(CTS)[sig][i]) # matplot
  plot_SS_Simulation(IC_new[[i]], simuBioTIP_g[[i]], 
                     main = paste("Delta Ic*"), 
                     ylab=NULL)
  
 }
dev.off()

  
# ######  old IC score, shulffing genes ##############
# set.seed(2021)
# IC_old <- simuIc_g <- list()  
# for(i in 1:length(sig)){
#   n <- length(CTS[sig][[i]])
#   IC_old[[i]] <- getIc(logmat, samplesL, CTS[sig][[i]], fun="cor", shrink=TRUE)
#   
#   simuIc_g[[i]]  <- simulation_Ic(n, samplesL, logmat, B=C,
#                                       fun="cor", shrink=TRUE)
# }
# names(simuIc_g) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
# names(IC_old) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
# 
# save(IC_old, simuIc_g,   file="simuIC_g.RData", compress=TRUE)
# 
# names(CTS)[sig]
# #[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2.5" "Day2.5"  
# 
# pdf(file="IC_vsSimulation_6CTS.pdf", height=10)
# par(mfrow=c(3,4))
# for(i in 1:length(IC_old)){
#   plot_Ic_Simulation(IC_old[[i]], simuIc_g[[i]], las = 2, ylab="old Ic",
#                      main= names(IC_new)[i],
#                      fun="boxplot", which2point= names(CTS)[sig][i]) # matplot
#   plot_SS_Simulation(IC_old[[i]], simuIc_g[[i]], 
#                      main = paste("Delta Ic"), 
#                      ylab=NULL)
#   
# }
# dev.off()

####### Verifying Tipping Point and using IC*  score  #################
######  BioTIP score, Shuffling label ##############
C= 1000
lengths(samplesL)
# Day0   Day1 Day1.5   Day2 Day2.5   Day3 
# 231    166     93    211    122    106 

set.seed(2021)
IC_new <- simuBioTIP_s <- list()  

for(i in 1:length(sig)){
  IC_new[[i]] <- getIc(logmat, samplesL, CTS[sig][[i]], fun="BioTIP", shrink=TRUE)
  simuBioTIP_s[[i]] <- matrix(nrow=length(samplesL), ncol=C)
 # rownames(simuBioTIP_s[[i]]) = names(CTS)[sig][i]
  rownames(simuBioTIP_s[[i]]) = names(samplesL)
  for(j in 1:length(samplesL)) { 
    #ns <- length(samplesL[names(CTS[sig])][[i]])  # for this state only
    ns <- length(samplesL[[j]])  # for each state rewpectively 
    simuBioTIP_s[[i]][j,]  <- simulation_Ic_sample(logmat, ns, Ic=IC_new[[i]],
                                             genes=CTS[sig][[i]], B=C,
                                             fun="BioTIP", shrink=TRUE)
  }
}
names(simuBioTIP_s) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
names(IC_new) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')

save(simuBioTIP_s,   file="simuBioTIP_s.RData", compress=TRUE)

names(CTS)[sig]
#[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2.5" "Day2.5"  
pdf(file="BioTIP_vs_Shufflinglabel_6CTS.pdf", height=10)
par(mfrow=c(3,4))
for(i in 1:length(IC_new)){
  plot_Ic_Simulation(IC_new[[i]], simuBioTIP_s[[i]], las = 2, ylab="Ic*",
                     main= names(IC_new)[i],
                     fun="boxplot", which2point= names(CTS)[sig][i]) # matplot
  plot_SS_Simulation(IC_new[[i]], simuBioTIP_s[[i]], 
                     main = paste("Delta Ic*"), 
                     ylab=NULL)
  
}
dev.off()


# ######  old IC score, shulffing labels ##############
# set.seed(2021)
# IC_old <- simuIc_s <- list() 
# 
# for(i in 1:length(sig)){
#   IC_old[[i]] <- getIc(logmat, samplesL, CTS[sig][[i]], fun="cor", shrink=TRUE)
#   simuIc_s[[i]] <- matrix(nrow=length(samplesL), ncol=C)
#   rownames(simuIc_s[[i]]) = names(samplesL)
#   for(j in 1:length(samplesL)) { 
#     ns <- length(samplesL[[j]])  # for each state rewpectively 
#     simuIc_s[[i]][j,]  <- simulation_Ic_sample(logmat, ns, Ic=IC_new[[i]],
#                                                  genes=CTS[sig][[i]], B=C,
#                                                  fun="cor", shrink=TRUE)
#  }
# }
# names(simuIc_s) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
# names(IC_old) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
# 
# save(IC_old, simuIc_s,   file="simuIc_s.RData", compress=TRUE)
# 
# names(CTS)[sig]
# #[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2.5" "Day2.5"  
# pdf(file="IC_vs_Shufflinglabel_6CTS.pdf", height=10)
# par(mfrow=c(3,4))
# for(i in 1:length(IC_old)){
#   plot_Ic_Simulation(IC_old[[i]], simuIc_s[[i]], las = 2, ylab="old Ic",
#                      main= names(IC_new)[i],
#                      fun="boxplot", which2point= names(CTS)[sig][i]) # matplot
#   plot_SS_Simulation(IC_old[[i]], simuIc_s[[i]], 
#                      main = paste("Delta Ic"), 
#                      ylab=NULL)
#   
# }
# dev.off()


