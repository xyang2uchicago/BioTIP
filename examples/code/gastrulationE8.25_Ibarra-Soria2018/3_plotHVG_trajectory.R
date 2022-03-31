
setwd('F:/projects/scRNA/results/IbarraSoria2018_MouseE8.25/E8.25.2018_robustness/')

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
load(file='../normalized.RData')
sce
# class: SingleCellExperiment 
# dim: 20483 19386 
# metadata(0):
# assays(3): counts normcounts logcounts
# rownames(20483): Xkr4 Gm1992 ... Csf2ra Gm21060
# rowData names(5): GENEID SYMBOL SEQNAME HVGs.10 HVGs.20
# colnames(19386): mixedMesoderm.a_1 mixedMesoderm.b_1 ... mixedMesoderm.b_715 blood_2079
# colData names(3): label celltype subcelltype
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

samplesL <- split(rownames(colData(sce)),f = colData(sce)$subcelltype)
length(samplesL)  # 33
pseudoOrder <- c('extraembryonicMesoderm', 'endothelial.a', 'endothelial.b', 'endothelial.c', 'endothelial.d','blood',
                 'mesodermProgenitors', 'presomiticMesoderm.b' , 'presomiticMesoderm.a', 'somiticMesoderm',
                 'mixedMesoderm.a', 'pharyngealMesoderm', 'mixedMesoderm.b', 'cardiac.a',  'cardiac.b' ,   'cardiac.c'
)

samplesL <- samplesL[pseudoOrder]
length(samplesL)  #[1] 16
sce <- sce[,unlist(samplesL)]
dim(sce) #20483 11039 

x <- apply(logcounts(sce.bk), 1, mean)
table(x>0)
# FALSE  TRUE 
#  1092 19391 
sce.bk <- sec.bk[which(x>0),]

x <- apply(logcounts(sce.bk), 1, mean)
hist(log(x), 100, main= '19391 genes with positive mean')
abline(v=-5,col='red')
dev.copy2pdf(file='cut.meanexpress.logcount.pdf')

table(log(x) > -5)
# FALSE  TRUE 
# 6688 12703 
sce <- sce[which(log(x) > -5), ]
dim(sce)
#[1] 12703 11039

# load global HVG --------------------------------------
sce.bk <- sce
load('sce_16subtype.RData')
dim(sce)
# 4000 11039 

HVG <- rownames(sce)
length(HVG) # 4000

sce <- sce.bk
rm(sce.bk)

# load local HVG -----------------------------
load('subcelltype/BioTIP_top0.1FDR0.05_optimized_local.HVG_selection.RData')
hvg <- lapply(testres, rownames) %>% unlist() %>% unique()
length(hvg)  #[1] 2294

## plot for feature selection ----------------------

x <- apply(logcounts(sce), 1, mean)
y <- apply(logcounts(sce), 1, sd)

mycol <- rep('light grey', nrow(sce))
mycol[which(rownames(sce) %in% HVG)] <- 'pink'
mycol[which(rownames(sce) %in% hvg)] <- 'red'
table(mycol)
#  light grey       pink        red 
#        8703       1706       2294 


pdf(file=paste0(length(hvg),
                ".variable_RTF.selection.pdf"))
par(mfrow=c(2,2))
plot(x,y, col=mycol, xlab='Average Expression in log', ylab='Standard Deviation', 
     pch=20) 
smoothScatter(x,y, col=mycol, xlab='Average Expression in log', ylab='Standard Deviation', 
              pch=20 , nbin = 100,
              # ,colramp = colorRampPalette(c("white", "grey95", 
              #                              "grey90", "grey70", "grey50", "grey40",
              #                              "grey30", "grey20", "grey10", "grey5"))
) 

dev.off()



################################
## plot trajectory
################################
# prepare customized colors -------------------
sce$subcelltype <- factor(as.vector(sce$subcelltype))
library(MouseGastrulationData)
head(EmbryoCelltypeColours)
length(EmbryoCelltypeColours) # 37

subcelltype <- levels(sce$subcelltype)
MyEmbryoCelltypeColours <- EmbryoCelltypeColours[1:length(subcelltype)]
names(MyEmbryoCelltypeColours) <- levels(sce$subcelltype)

# buildtrajectory ---------------------------
library(scater)
by.cluster <- aggregateAcrossCells(sce, ids=colData(sce)$subcelltype)
centroids <- reducedDim(by.cluster, "PCA")
colData(dimred.sce)

dmat <- dist(centroids)
dmat <- as.matrix(dmat)
g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
mst <- igraph::minimum.spanning.tree(g)


pdf(file="trajectory_UMAP.pdf", height=5)

#set.seed(1000)
plot(mst)

pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
coords <- reducedDim(by.cluster, "UMAP")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)

#gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="UMAP", colour_by='subcelltype', text_by='subcelltype',
                 text_colour=MyEmbryoCelltypeColours, 
                 point_size=0.2, text_siz=3) + xlim(-10,5) +
    scale_color_manual(values=MyEmbryoCelltypeColours)#,
  
  plotReducedDim(sce, dimred="UMAP", colour_by='subcelltype', text_by='subcelltype',
                 text_colour=MyEmbryoCelltypeColours, 
                 point_size=0.2, text_siz=3) + xlim(-10,5) +
    scale_color_manual(values=MyEmbryoCelltypeColours) + 
    geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))#,
#   ncol=2
# )

dev.off()



