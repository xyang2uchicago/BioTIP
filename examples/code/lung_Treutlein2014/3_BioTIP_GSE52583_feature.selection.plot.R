setwd('F:/projects/BioTIP/result/GSE52583/')

# This code does:
# (following the code of BioTIP_E8.25_mesoderm_cluster_PCCsNoShrink.R)
# plot RFT-selected genes: 29k (not shown) -> 10.9k expressed (grey) ->4k HGV (blue) -> 3073 RTF selected (red)

#-----------------------------------
# normalization

# normalized data were downloaded and reanalyzed
# refer to GSE87038_E8.25_normalize.R
#-----------------------------------------
#this document uses library size to continue on later analysis
#------------------------------------------

library(scater)
packageVersion('scater') # '1.18.3'
library(scran)
packageVersion('scran') # '1.18.1'
# 
# 
library(BioTIP)
packageVersion("BioTIP")  #[1] '1.4.0'

# 
###################################################################
## 22.9k measured transcripts -> 10.3k expressed genes -> 3.2k HGV
###################################################################
require(stringr)
require(psych)
require(igraph)
 
load(file="GSE52583_annot.cds_Clustered_thresholded.RData")
annot.cds
# CellDataSet (storageMode: environment)
# assayData: 22854 features, 196 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: GSM1271862 GSM1271863 ... GSM1272062 (196 total)
# varLabels: platform_id instrument_model ... State (21 total)
# varMetadata: labelDescription
# featureData
# featureNames: 0610005C13Rik 0610007C21Rik ... l7Rn6 (22854 total)
# fvarLabels: ensembl_gene_id gene_short_name ... use_for_ordering (11
#                                                                   total)
# fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: 
# 
table(fData(annot.cds)$use_for_ordering)
# FALSE  TRUE 
# 12495 10359


annot.cds <- annot.cds[which(fData(annot.cds)$use_for_ordering),]
dim(annot.cds)
#Features  Samples 
# 10359      196

# select coding genes, lncRNAs, and miRNAs for the downstream analysis ------------------------
table(fData(annot.cds)$use_for_ordering, fData(annot.cds)$gene_biotype %in% c('lncRNA', 'miRNA','protein_coding') )
#       FALSE  TRUE
# FALSE   380 12115
# TRUE    108 10251

x <- which(fData(annot.cds)$use_for_ordering &
             fData(annot.cds)$gene_biotype %in% c('lncRNA', 'miRNA','protein_coding') )
length(x) # 10251
vars <- rowVars(log2(exprs(annot.cds[x,])+1))
hist(vars, 100)
table(vars>0.5)
# FALSE  TRUE 
#  7053  3198 
HVG <- rownames(annot.cds)[x][which(vars>0.5)]



#########################################################################
## RTF selected local HVGs for 131 cells along the AT2 trajectory
## demonstrate the BioTIP results using the Leiden_0.4 clusters 
######################################################################## 
load('BioTIP_GSE52583_robustness/AT2.sce.RData') 
sce
# class: SingleCellExperiment 
# dim: 3198 131 
# metadata(0):
#   assays(1): logFPKM
# rownames(3198): 0610007C21Rik 0610007L01Rik ... Zyx l7Rn6
# rowData names(11): ensembl_gene_id gene_short_name ...
# num_cells_expressed use_for_ordering
# colnames(131): E18_1_C08 E18_1_C09 ... E16_1_C94 E16_1_C95
# colData names(17): age cellName ... C_SNNGraph_k8 C_Soft
# reducedDimNames(2): PCA TSNE
# altExpNames(0):
  
fid <- "BioTIP_GSE52583_robustness/C_Leiden_0.4/"
load(file=paste(fid,"BioTIP_top0.1FDR0.2_optimized_local.HVG_selection.RData"))
names(testres)
#[1] "1" "2" "3"

table(colData(sce)$C_Leiden_0.4)
# 1  2  3  4  5 
# 46 32 24 17 12 


########################
## plot
#############################
dat <- log2(exprs(annot.cds[,colData(sce)$GEO.access])+1)
dim(dat)
#[1] 10359   131

x <- apply(dat, 1, mean)
y <- apply(dat, 1, sd)

# mycol <- rep('light grey', nrow(dec.pois))
# mycol[which(rownames(dec.pois) %in% rownames(sce))] <- 'grey0'
mycol <- rep('light grey', nrow(dat))
mycol[which(rownames(dat) %in% HVG)] <- 'pink'
mycol[which(rownames(dat) %in% unique(unlist(lapply(testres, rownames))))] <- 'red'
(z <- table(mycol))
# light grey       pink        red 
#       7161       2444        754


pdf(file=paste0(fid, length(HVG),
                ".variable_RTF.selection.pdf"))
par(mfrow=c(2,2))
plot(x,y, col=mycol, xlab='Average Expression in log2', ylab='Standard Deviation', 
     pch=20) 
smoothScatter(x,y, col=mycol, xlab='Average Expression in log2', ylab='Standard Deviation', 
              pch=20 #, nbin = 100,
              # ,colramp = colorRampPalette(c("white", "grey95", 
              #                              "grey90", "grey70", "grey50", "grey40",
              #                              "grey30", "grey20", "grey10", "grey5"))
              ) 

# plot(dec.pois$mean, dec.pois$bio, col=mycol, xlab='Average Expression in log2', 
#      ylab='Biological Variance', 
#      pch=20) 
# suppressWarnings(
#   smoothScatter(dec.pois$mean, dec.pois$bio, col=mycol, xlab='Average Expression in log2', 
#                 ylab='Biological Variance', 
#      pch=20 #,nbin = 100
#      ) 
# )
dev.off()

