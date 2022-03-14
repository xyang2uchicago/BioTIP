## NOTE that some of the intermediately resultant R objects were not pushed to github yet and are available on request.
## This code only tested on local working fold yet.

setwd("F:/projects/scRNA/results/AJ/GSE130146_xy/Results_1k/LibrarySize/GSE130146_robustness/") 

# PubMed ID	29311656
# GRCm38.p4
# annotation from the Ensembl database, version 84.
#
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)


###########################################################
## section 1) plot feature selection:
## 20483 measured genes
## 12703 expressed genes (average logcounts > 10^(-5))
## 4000 global HVG
## 2294 local hvg selected by BioTIP
###########################################################

# load sce-------------------
load(file="../dimred.sce_afterPCA.RData")
dimred.sce 
# class: SingleCellExperiment 
# dim: 33456 1731 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(33456): ENSMUSG00000064842 ENSMUSG00000051951 ... ENSMUSG00000096730 ENSMUSG00000095742
# rowData names(2): ID Symbol
# colnames(1731): AAACCTGAGTTTCCTT-1 AAACCTGCATGAAGTA-1 ... TTTGTCATCGGCTTGG-1 TTTGTCATCTGCTTGC-1
# colData names(4): Sample Barcode sizeFactor label
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):
  
#dropping endoderm
dropcols = colData(dimred.sce)$label == 16 |colData(dimred.sce)$label == 2  |
  colData(dimred.sce)$label == 15 |colData(dimred.sce)$label == 13
dimred.sce = dimred.sce[,!dropcols]

dim(dimred.sce)
#[1] 33456  1531

x <- apply(logcounts(dimred.sce), 1, mean)
table(x>0)
# FALSE  TRUE 
# 15411 18045
sce.bk <- dimred.sce[which(x>0),]
rownames(sce.bk) <- rowData(sce.bk)$Symbol

x <- apply(logcounts(sce.bk), 1, mean)
hist(log(x), 100, main= '15411 genes with positive mean')
abline(v=-6,col='red')
dev.copy2pdf(file='cut.meanexpress.logcount.pdf')

table(log(x) > -6)
# FALSE  TRUE 
# 2845 15200 
sce.bk <- sce.bk[which(log(x) > -6), ]
dim(sce.bk)
#[1] 15200   1531

# load global HVG --------------------------------------
load('sce.GSE130146_noenderdormPgerm.RData')
dim(sce)
# 4000 1531

HVG <- rownames(sce)
length(HVG) # 4000

sce <- sce.bk
rm(sce.bk)

# load local HVG -----------------------------
load('C_SNNGraph_k5/BioTIP_top0.1FDR0.2_optimized_local.HVG_selection.RData')
hvg <- lapply(testres, rownames) %>% unlist() %>% unique()
length(hvg)  #[1] 2078

## plot for feature selection ----------------------

x <- apply(logcounts(sce), 1, mean)
y <- apply(logcounts(sce), 1, sd)

mycol <- rep('light grey', nrow(sce))
mycol[which(rownames(sce) %in% HVG)] <- 'pink'
mycol[which(rownames(sce) %in% hvg)] <- 'red'
table(mycol)
# light grey       pink        red 
#     13357        882        961 

pdf(file=paste0(length(hvg),
                ".variable_RTF.selection.pdf"))
par(mfrow=c(2,2))
plot(x,y, col=mycol, xlab='Average Expression in log', ylab='Standard Deviation', 
     pch=20) 
smoothScatter(x,y, col=mycol, xlab='Average Expression in log', ylab='Standard Deviation', 
              pch=20 , nbin = 100
) 

dev.off()



################################
## plot trajectory
################################

## Minimum spanning tree constructed trajectory ---------------
library(scater)
by.cluster <- aggregateAcrossCells(dimred.sce, ids=colData(dimred.sce)$label)
centroids <- reducedDim(by.cluster, "PCA")

dmat <- dist(centroids)
dmat <- as.matrix(dmat)
set.seed(1000)
g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
mst <- igraph::minimum.spanning.tree(g)

pdf(file="EB_trajectory_libsf_without_endodermPgerm.pdf")  #!!!!!!!!!!!!!!!!!!!!
plot(mst)

pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
coords <- reducedDim(by.cluster, "TSNE")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)

plotTSNE(dimred.sce, colour_by="label", 
         text_by="label", text_size = 8)  

plotTSNE(dimred.sce, colour_by="label", 
         text_by="label", text_size = 8) + 
  geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))

dev.off()



####################################################################
## plot the observed Ic.shrink vs random Ic.shrink in a better view
####################################################################
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/myplotIc.R')

setwd("F:/projects/scRNA/results/AJ/GSE130146_xy/Results_1k/LibrarySize/GSE130146_robustness/C_SNNGraph_k5")

load('BioTIP.res.RData')
CTS.candidate = res$CTS.candidate 
BioTIP_scores = res$Ic.shrink 
load('BioTIP_top0.1FDR0.2_IC_sim.Permutation.RData')


n.permutation=100
local.IC.p=TRUE
n <- length(CTS.candidate)
p.IC <- rep(1, n)
if(length(SimResults_g)>0) {
  SimResults_4_p.IC = SimResults_g } else if(length(SimResults_s)>0) {
    SimResults_4_p.IC = SimResults_s 
  } else if(length(SimResults_b)>0) {
    SimResults_4_p.IC = SimResults_b
  }
if(length(SimResults_4_p.IC)>0) {
  for(i in 1:n){
    interesting =  names(BioTIP_scores[i]) ## uodated 1/30/2022
    #first p value calculated for exactly at tipping point
    p = length(which(SimResults_4_p.IC[[i]][interesting,] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
    p = p/n.permutation  # ncol(SimResults_4_p.IC[[i]])
    #second p value calculated across all statuses
    p2 = length(which(SimResults_4_p.IC[[i]] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
    p2 = p2/n.permutation  #ncol(SimResults_4_p.IC[[i]])
    p.IC[i] = p2/nrow(SimResults_4_p.IC[[i]])
  }
}   

#outputpath = paste0(BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
#pdf.file.name = paste0(outputpath,"IC_Delta_SimresultGene.pdf")
pdf.file.name = 'BioTIP_top0.1FDR0.2_IC_Delta_SimresultGene.pdf'
myplotIc(pdf.file.name, BioTIP_scores, CTS.candidate, SimResults_g, width=12, height=4.5, nn=6 )
