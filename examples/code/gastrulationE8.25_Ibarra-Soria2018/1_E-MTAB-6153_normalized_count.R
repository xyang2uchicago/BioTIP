
setwd('F:/projects/scRNA/results/IbarraSoria2018_MouseE8.25/')

# PubMed ID	29311656
# GRCm38.p4
# annotation from the Ensembl database, version 84.
#
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(scater)
library(scran)
library(pheatmap)
library(limma)
library(edgeR)
library(dynamicTreeCut)
library(cluster)

library(MouseGastrulationData)
head(EmbryoCelltypeColours)
length(EmbryoCelltypeColours) # 37

##################################################
## section 1) load rwas and normalized data
##################################################

## http://bioinformatics.age.mpg.de/presentations-tutorials/presentations/modules/single-cell//bioconductor_tutorial.html
countsfile="../../data/IbarraSoria2018_MouseE8.25/E-MTAB-6153.processed.2/rawCounts.tsv"
all.counts <- as.data.frame(read.delim(countsfile, header=T, row.names=1))
all.counts[1:3,1:5]
# mixedMesoderm.a_1 mixedMesoderm.b_1 neuralCrest_1
# ENSMUSG00000051951                 0                 0             0
# ENSMUSG00000089699                 0                 0             0
# ENSMUSG00000025900                 0                 0             0
# extraembryonicEndoderm_1 badQuality_1
# ENSMUSG00000051951                        0            0
# ENSMUSG00000089699                        0            0
# ENSMUSG00000025900                        0            0
dim(all.counts)
#[1] 20809 20819
# rownames(all.counts) <- all.counts$ID
# all.counts <- as.matrix(all.counts[,-1])
#For convenience, the counts for spike-in transcripts and endogenous genes are stored in a SingleCellExperiment object from the SingleCellExperiment package.

# droplet based single-cell RNA-sequencing to address this by profiling ~20000 cells from C57BL/6 E8.25 mouse embryos.
# What has been normalized by the authors were:
# The data were normalized for cell-specific biases using a
# method previously proposed33 and implemented in the Bioconductor package
# scran34. To calculate size factors, genes with a mean expression < 0.1 were filtered
# out; the quickCluster function was used to obtain the initial clustering of the
# cells (method igraph). The estimated size factors were used to normalize all
# genes expressed in at least one cell. Normalized counts are provided with the
# ArrayExpress submission.

dat <- read.csv('../../data/IbarraSoria2018_MouseE8.25/normalisedCounts.tsv', row.name=1, head=T, sep='\t')
dim(dat)  # 20806 19386
head(dat)    # [1] 20,806 19,386
class(dat)  #[1] "data.frame"


# all assays must have the same nrow and ncol
all(rownames(dat) %in% rownames(all.counts)) T
all(colnames(dat) %in% colnames(all.counts)) T

sce <- SingleCellExperiment(list(counts=all.counts[rownames(dat), colnames(dat)], normcounts=dat))
dim(sce) # 20806 19386
rm(all.counts, dat)



library(AnnotationHub)
# the author used Gencode GRCh38.primary (ftp://ftp.sanger.ac.uk/pub/
# gencode/Gencode_human/release_25/GRCh38.primary_assembly.genome.fa.gz)
# as the genome sequence and Gencode v25 transcript annotations 
tmp <- AnnotationHub()
#snapshot date: 2020-10-27

query(tmp, pattern = c("Mus musculus", "EnsDb"))  # Homo sapiens
## snapshotDate(): 2020-10-27
# # $dataprovider: Ensembl
# # $species: Mus musculus
# # $rdataclass: EnsDb
# # additional mcols(): taxonomyid, genome, description,
# #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
# #   rdatapath, sourceurl, sourcetype 
# # retrieve records with, e.g., 'object[["AH53222"]]' 
# 
# title                             
# AH53222 | Ensembl 87 EnsDb for Mus Musculus 
# AH53726 | Ensembl 88 EnsDb for Mus Musculus 
# AH56691 | Ensembl 89 EnsDb for Mus Musculus 
# AH57770 | Ensembl 90 EnsDb for Mus Musculus 
# AH60788 | Ensembl 91 EnsDb for Mus Musculus 
# ...       ...                               
# AH78811 | Ensembl 99 EnsDb for Mus musculus 
# AH79718 | Ensembl 100 EnsDb for Mus musculus
# AH83247 | Ensembl 101 EnsDb for Mus musculus
# AH89211 | Ensembl 102 EnsDb for Mus musculus
# AH89457 | Ensembl 103 EnsDb for Mus musculus
ens.mm.v103 <- AnnotationHub()[["AH89457"]]
class(ens.mm.v103)
#"EnsDb"  
anno <- select(ens.mm.v103, keys=rownames(sce), 
               keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))

rowData(sce) <- anno[match(rownames(sce), anno$GENEID),]

chr.loc <- mapIds(ens.mm.v103, keys=rownames(sce),
                  keytype="GENEID", column="SEQNAME")
# Unable to map 309 of 20806 requested IDs. 
is.mito <- which(chr.loc=="MT")
# is.mito <- (chr.loc=="MT")
# summary(is.mito)
# # Mode   FALSE    TRUE    NA's 
# # logical   20489      11     309 

is.spike <- grepl("^ERCC", rownames(sce))
sum(is.spike) # 0



########################################
## quality control
########################################

library(scater)
df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
df
#DataFrame with 1444 rows and 6 columns

# for filtering the cells, we refer to the downloaded normalized data matrix.
# Here, we verify what the paper said that 
# We removed all cells that expressed < 1,000 genes or that had > 3%
# of their transcripts mapped to mitochondrial genes. We further removed any cells
# that expressed both Xist and any of Kdm5d, Eif2s3y, Gm29650, Uty or Ddx3y (genes
# in the Y chromosome), as these are probably doublets. We identified 400 cells that
# could be affected by index swapping (as they share the same cell barcode with
# another cell), even though the rates of this phenomenon are very low for the HiSeq 2500.
qc.nexprs <- df$detected < 1000
qc.mito <- df$subsets_Mito_percent > 3
# tmp <- rbind(counts(sce)[which(rowData(sce)$SYMBOL == 'Xist'),], 
#              apply(counts(sce)[which(rowData(sce)$SYMBOL %in% c('Kdm5d', 'Eif2s3y', 'Gm29650', 'Uty', 'Ddx3y')),],2,sum)>1 )
# qc.chrY <- apply(tmp,2,sum)>1
discard <- qc.nexprs | qc.mito  #| qc.chrY #| qc.chr
DataFrame(NExprs=sum(qc.nexprs), MitoProp=sum(qc.mito), # qc.chrY=sum(qc.chrY), Loc=sum(qc.chr),
         Total=sum(discard))
#        NExprs  MitoProp     Total
#      <integer> <integer> <integer>
# 1         0         0         0

# filted <- sce[,!discard]
# dim(filted)


# removing doublets -------------------
# CellRanger version 3 automatically performs cell calling using an algorithm similar to emptyDrops().
## However, this data was generated using CellRanger v2.0.0.
library(DropletUtils)
bcrank <- barcodeRanks(counts(sce))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

metadata(bcrank)$inflection #[1] 3848
metadata(bcrank)$knee #[1] 11242
dev.copy2pdf(file='UMO_count_each_barcode.pdf')


# emptyDrops performs Monte Carlo simulations to compute p-values,
# so we need to set the seed to obtain reproducible results.
set.seed(100)
limit <- metadata(bcrank)$inflection
e.out <- emptyDrops(counts(sce), lower=limit)

#Ideally, the distribution should be close to uniform -- which don't hold true for this case!!
# Large peaks near zero indicate that barcodes with total counts below lower are not all ambient in origin.
hist(e.out$PValue, 100) 

#The Limited field in the output indicates whether or not the computed  
#p-value for a particular barcode is bounded by the number of iterations. 
#If any non-significant barcodes are TRUE for Limited, we may need to increase the number of iterations. 
table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)
#         Limited
# Sig     FALSE TRUE
#  FALSE   969     0
#  TRUE    547 17852


set.seed(100)
#assumes that barcodes with low total UMI counts are empty droplets. 
# Thus, the null hypothesis should be true for all of these barcodes. 
# We can check whether the hypothesis testing procedure holds its size 
# by examining the distribution of p-values for low-total barcodes with (test.ambient=TRUE). Ideally, the distribution should be close to uniform (Figure 15.2)
all.out <- emptyDrops(counts(sce), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0], 100,
     xlab="P-value", main="", col="grey80") 
dev.copy2pdf(file=paste0('hist_assumed_empty_droplets.pdf'))

table(e.out$FDR <= 0.001) # not appliable, eitherwise 18k out of 19k genes will be excluded !!!
# FALSE  TRUE 
#  969 18399 

# filted <- sce[,which(e.out$FDR <= 0.001 & !discard)]
# dim(filted)



# Removing ambient contamination (generally be omitted from most analyses)
# The removeAmbience() function from DropletUtils will remove the contamination from the cluster-level profiles 
# and propagate the effect of those changes back to the individual cells.

# Removing swapped molecules(not appliable for Hiseq 500)
# Some of the more recent DNA sequencing machines released by Illumina 
# (e.g., HiSeq 3000/4000/X, X-Ten, and NovaSeq) use patterned flow cells 
# to improve throughput and cost efficiency. However, in multiplexed pools, 
# the use of these flow cells can lead to the mislabelling of DNA molecules 
# with the incorrect library barcode (Sinha et al. 2017), 
# a phenomenon known as "barcode swapping".
#
# library(DropletTestFiles)
# swap.files <- listTestFiles(dataset="bach-mammary-swapping")
# swap.files <- swap.files[dirname(swap.files$file.name)=="hiseq_4000",]
# swap.files <- vapply(swap.files$rdatapath, getTestFile, prefix=FALSE, "")
# names(swap.files) <- sub(".*_(.*)\\.h5", "\\1", names(swap.files))


# for filter genes
# removing those mapped to abnormal chromosome
table(chr.loc %in% c(1:22,'X','Y'))
# FALSE  TRUE 
# 323 20483 
sce <- sce[which(chr.loc %in% c(1:22,'X','Y')),]
dim(sce)
#  20483 19386


## Doublet detection
# doublets are artifactual libraries generated from two cells. 
# They typically arise due to errors in cell sorting or capture, 
# especially in droplet-based protocols (Zheng et al. 2017) involving thousands of cells.
# Like 'findMarkers', this function will automatically
# retrieve cluster assignments from 'colLabels'.
# the parameter clusters taken from colLabels(x) by default; need at least three clusters to detect doublet clusters
#library(scDblFinder)


# library(scran)
# markers <- findMarkers(sce, direction="up")
# dbl.markers <- markers[[chosen.doublet]]
# 
# library(scater)
# chosen <- rownames(dbl.markers)[dbl.markers$Top <= 10]
# plotHeatmap(sce.mam, order_columns_by="label", features=chosen, 
#             center=TRUE, symmetric=TRUE, zlim=c(-5, 5))

assayNames(sce)
# [1] "counts"     "normcounts"

########################################
## log transform the normalized counts

logcounts(sce) <- log2(normcounts(sce) + 1)   # by default was 1
assayNames(sce)
#[1] "counts"     "normcounts" "logcounts" 

# save(sce,file='filtered.RData')



#########################
## feature selection
#########################
library(scran)
dec <- modelGeneVar(sce)

# fit a trend to the variance with respect to abundance across all genes
# Visualizing the fit:
fit <- metadata(dec)
plot(fit$mean, fit$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)
dev.copy2pdf(file='fit_variance.pdf')

rownames(sce) <- uniquifyFeatureNames(rowData(sce)$GENEID, rowData(sce)$SYMBOL)

# Ordering by most interesting genes for inspection.
dec[order(dec$bio, decreasing=TRUE),] 

chosen <- getTopHVGs(dec, prop=0.1)
str(chosen)
length(chosen)
#[1]  1355
chosen[1:4]
# [1] "ENSMUSG00000052217" "ENSMUSG00000055609" "ENSMUSG00000061808"
# [4] "ENSMUSG00000069919"

rowSubset(sce,'HVGs.10') <- rowData(sce)$GENEID %in% chosen # stored in the default 'subset'.
rowSubset(sce, "HVGs.20") <- rowData(sce)$GENEID %in% getTopHVGs(dec, prop=0.2)


################################
## Dimensionality reduction, based on the top 10% HGVs
################################
set.seed(100) # See below.
sce <- runPCA(sce, subset_row=rowData(sce)$HVGs.10) 

## chose the number of PCs
percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")

# Percentage of variance explained is tucked away in the attributes.
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow  # 5
plot(percent.var, xlab="PC", ylab="Variance explained (%)")
abline(v=chosen.elbow, col="red")


# ## using the technical noise
# library(scran)
# set.seed(111001001)
# denoised <- denoisePCA(sce, technical=dec, subset.row=rowData(sce)$HVGs.10)
# ncol(reducedDim(denoised))  # 5


## using randomized matrix
set.seed(100010)
horn <- PCAtools::parallelPCA(logcounts(sce)[which(rowData(sce)$HVGs.10),],
                              BSPARAM=BiocSingular::IrlbaParam(), niters=10)
horn$n
## [1] 21

pcs <- reducedDim(sce)
choices <- getClusteredPCs(pcs)
val <- metadata(choices)$chosen

pdf('chose_number_of_PCA.pdf')

plot(percent.var, xlab="PC", ylab="Variance explained (%)", main=outputs[i])
abline(v=chosen.elbow, col="red")

plot(horn$original$variance, type="b", log="y", pch=16)
permuted <- horn$permuted
for (j in seq_len(ncol(permuted))) {
  points(permuted[,j], col="grey80", pch=16)
  lines(permuted[,j], col="grey80", pch=16)
}
abline(v=horn$n, col="red")

#  The red unbroken line represents the theoretical upper constraint on the number of clusters, 
# while the grey dashed line is the number of PCs suggested by getClusteredPCs()
plot(choices$n.pcs, choices$n.clusters,
     xlab="Number of PCs", ylab="Number of clusters", main=outputs[i])
abline(a=1, b=1, col="red")
abline(v=val, col="grey80", lty=2)

dev.off()

####################
## using the top 21 PCS !!!!!!!!!!!!!
dim(reducedDim(sce, "PCA")) 
#  19386    50
reducedDim(sce, "PCA") <- reducedDim(sce, "PCA")[,1:horn$n]



set.seed(00101001101)

# runTSNE() stores the t-SNE coordinates in the reducedDims
# for re-use across multiple plotReducedDim() calls.
# Note that by default, scater::runTSN using the parameter n.componets=50
# However, only the 21 PCs were recorded !!!!!!!!!!!
sce <- runTSNE(sce, dimred="PCA")
#plotReducedDim(sce, dimred="TSNE", colour_by=clust)

sce <- runUMAP(sce, dimred="PCA" )
#plotReducedDim(sce, dimred="UMAP", colour_by=clust)


###############################
## cluster cells, buildSNNGraph, type = "jaccard"
# the paper said that
# cells with similar transcriptional profiles were clustered into
# 33 different groups, as indicated by the different colours. 
# Each cluster was annotated based on the expression of marker genes in 20 different major cell types.
# Several cell types are composed of two or more clusters.
##############################
library(scran)
# type="number" will weight edges based on the number of nearest neighbors that are shared between two cells. 
# type="jaccard" will weight edges according to the Jaccard index of the two sets of neighbors.
g <- buildSNNGraph(sce, k=20, use.dimred = 'PCA', type="jaccard")
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
#clust
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 577  365 1008 1034  884 1060  818 2329  411  372 1114  911 1252 1157  613 1847 
# 17   18   19   20   21   22   23   24   25   26   27   28   29 
# 222  593  263  741  199  273  152  253  205  307  322   80   24


colLabels(sce) <- factor(clust)

library('cluster')
library(bluster)
# dist <- dist(reducedDim(sce, "PCA"))
# sil <- silhouette(clust, dist = dist)
# plot(sil)

MyEmbryoCelltypeColours <- EmbryoCelltypeColours[1:length(unique(clust))]
names(MyEmbryoCelltypeColours) <- 1:length(unique(clust))

pdf(file='ReduceDim_SNNGraph_20_jaccard.pdf', height=4.5, width=10) #!!!!!!!!!!!!!!!!!

# Performing the calculations on the PC coordinates, like before.
sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=clust)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, clust, sil.data$other))
sil.data$cluster <- factor(clust)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="PCA", colour_by='label', text_by='label',point_size=0.2) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "TSNE"), clusters=clust)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, clust, sil.data$other))
sil.data$cluster <- factor(clust)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="TSNE", colour_by='label', text_by='label',point_size=0.2) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "UMAP"), clusters=clust)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, clust, sil.data$other))
sil.data$cluster <- factor(clust)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="UMAP", colour_by='label', text_by='label',point_size=0.2) +
     scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ncol=2
)

dev.off()

# sce.iPS[[i]] <- sce 
# save(sce, file='normalized.RData')


# celltypePlot <- function(sce, dirmed='UMAP',pch = 19, col2plot='label')
# {
#   # plot(
#   #   x = reducedDim(sce, type=dirmed)[, 1],
#   #   y = reducedDim(sce, type=dirmed)[, 2],
#   #   col = EmbryoCelltypeColours[colData(sce)[,col2plot]],
#   #   pch = pch, cex=0.5,
#   #   xaxt = "n", yaxt = "n",
#   #   xlab = paste0(dirmed,'1'), ylab = paste0(dirmed,'2')
#   # )
#      x = reducedDim(sce, type=dirmed)[, 1]
#      y = reducedDim(sce, type=dirmed)[, 2]
#      ggplot(data.frame(reducedDim(sce,dirmed)[,1:2]),
#           aes(x=X1,y=X2, color = EmbryoCelltypeColours)) +
#      geom_point()
# }


####################################################################
## annotate the cluster used 
## the paper said that
## Pearson's chi-squared test (P values corrected for multiple
## testing using the Benjamini and Hochberg method; Supplementary Table 1
####################################################################

## extract the author's annotations
tmp <- colnames(sce)
subcelltype <- unlist(lapply(tmp, function(x) unlist(strsplit(x, '_'))[1]))
length(unique(subcelltype))  # 33
celltype <- unlist(lapply(subcelltype, function(x) unlist(strsplit(x, '.', fixed=T))[1]))
length(unique(celltype))  # 20

df <- table(sce$label, celltype)
as.data.frame(df)
#      Var1               celltype Freq
# 1      1                 amnion    0
# 2      2                 amnion    0
# 3      3                 amnion    4
# 4      4                 amnion    0
# 5      5                 amnion   16
# 6      6                 amnion    0
df2 <- table(sce$label, subcelltype)

gridExtra::grid.arrange(
  ggplot(as.data.frame(df), aes(x=Var1, y=Freq, fill=as.factor(celltype) )) + 
    geom_bar(stat = "identity"),
  ggplot(as.data.frame(df2), aes(x=Var1, y=Freq, fill=as.factor(subcelltype) )) + 
    geom_bar(stat = "identity"),
  ncol=1
)
dev.copy2pdf(file='bar_cluster_vs_celltype.pdf')


pdf(file='bar_celltype_vs_cluster.pdf', height=10, width=10)
gridExtra::grid.arrange(
  ggplot(as.data.frame(df), aes(x=celltype, y=Freq, fill=as.factor(Var1) )) + 
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 90)),
  ggplot(as.data.frame(df2), aes(x=subcelltype, y=Freq, fill=as.factor(Var1) )) + 
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90)),
  ncol=1
)
dev.off()


sce$celltype <- as.factor(celltype)
sce$subcelltype <- as.factor(subcelltype)
save(sce, file='normalized.RData')  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 
# plotExpression(sce, features=c("ETV2", "TAL1", 'SOX9','KDR'),
#                x=I(sce$label))
# 
# dev.copy2pdf(file='expression_marker.pdf')
# 

MyEmbryoCelltypeColours <- EmbryoCelltypeColours[1:length(unique(celltype))]
names(MyEmbryoCelltypeColours) <- unique(celltype)

pdf(file='ReduceDim_celltype.pdf', height=4.5, width=10) #!!!!!!!!!!!!!!!!!

# Performing the calculations on the PC coordinates, like before.
sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=celltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, celltype, sil.data$other))
sil.data$cluster <- factor(celltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="PCA", colour_by='celltype', text_by='celltype', 
                 point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90)),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "TSNE"), clusters=celltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, celltype, sil.data$other))
sil.data$cluster <- factor(celltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="TSNE", colour_by='celltype', text_by='celltype',
                 point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90)),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "UMAP"), clusters=celltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, celltype, sil.data$other))
sil.data$cluster <- factor(celltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="UMAP", colour_by='celltype', text_by='celltype', 
                 point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90)),
  ncol=2
)

dev.off()


MyEmbryoCelltypeColours <- EmbryoCelltypeColours[1:length(unique(subcelltype))]
names(MyEmbryoCelltypeColours) <- unique(subcelltype)

pdf(file='ReduceDim_subcelltype.pdf', height=4.5, width=18) #!!!!!!!!!!full dataset, not used !!!!!!!

# Performing the calculations on the PC coordinates, like before.
sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=subcelltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, subcelltype, sil.data$other))
sil.data$cluster <- factor(subcelltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="PCA", colour_by='subcelltype', text_by='subcelltype', 
                 point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95)),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "TSNE"), clusters=subcelltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, subcelltype, sil.data$other))
sil.data$cluster <- factor(subcelltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="TSNE", colour_by='subcelltype', text_by='subcelltype',
                 point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95)),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "UMAP"), clusters=subcelltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, subcelltype, sil.data$other))
sil.data$cluster <- factor(subcelltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="UMAP", colour_by='subcelltype', text_by='subcelltype', 
                 point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95)),
  ncol=2
)

dev.off()

#########################################################
# UMAP of 16 subcelltypes of developing mesoderm cells, USED !!!!!!!!!!!!!!!
################################################
load(file='F:/projects/scRNA/results/IbarraSoria2018_MouseE8.25/normalized.RData')
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
lengths(samplesL)
# extraembryonicMesoderm          endothelial.a          endothelial.b          endothelial.c
# 1017                    239                    128                    458
# endothelial.d                  blood    mesodermProgenitors   presomiticMesoderm.b
# 47                   2079                   1696                    275
# presomiticMesoderm.a        somiticMesoderm        mixedMesoderm.a     pharyngealMesoderm
# 875                    739                   1328                    743
# mixedMesoderm.b              cardiac.a              cardiac.b              cardiac.c
# 715                    234                    184                    282

sce <- sce[,unlist(samplesL)]
sce
# dim: 20483 11039
sce$celltype <- as.factor(as.vector(sce$celltype))
sce$subcelltype <- as.factor(as.vector(sce$subcelltype))
sce$label <- as.factor(as.vector(sce$label))

MyEmbryoCelltypeColours <- EmbryoCelltypeColours[1:length(unique(sce$subcelltype))]
names(MyEmbryoCelltypeColours) <- levels(sce$subcelltype)


library(bluster)
subcelltype = sce$subcelltype


pdf(file='ReduceDim_subcelltype_16mesoderm.pdf', height=4.5, width=18) #!!!!!!!!!!full dataset, not used !!!!!!!

# Performing the calculations on the PC coordinates, like before.
sil.approx <- approxSilhouette(reducedDim(sce, "PCA"), clusters=subcelltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, subcelltype, sil.data$other))
sil.data$cluster <- factor(subcelltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="PCA", colour_by='subcelltype', text_by='subcelltype', 
                 text_colour=MyEmbryoCelltypeColours, point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    #scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95)),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "TSNE"), clusters=subcelltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, subcelltype, sil.data$other))
sil.data$cluster <- factor(subcelltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="TSNE", colour_by='subcelltype', text_by='subcelltype',
                 text_colour=MyEmbryoCelltypeColours, point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    #scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95)),
  ncol=2
)

sil.approx <- approxSilhouette(reducedDim(sce, "UMAP"), clusters=subcelltype)
sil.data <- as.data.frame(sil.approx)
sil.data$closest <- factor(ifelse(sil.data$width > 0, subcelltype, sil.data$other))
sil.data$cluster <- factor(subcelltype)

gridExtra::grid.arrange(
  plotReducedDim(sce, dimred="UMAP", colour_by='subcelltype', text_by='subcelltype', 
                 text_colour=MyEmbryoCelltypeColours, point_size=0.2, text_siz=3) +
    scale_color_manual(values=MyEmbryoCelltypeColours),
  ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) + 
    ggbeeswarm::geom_quasirandom(method="smiley")+
    #scale_color_manual(values=MyEmbryoCelltypeColours)+
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95)),
  ncol=2
)

dev.off()



#############################3
## biomarkers, for the 16 mesoderm celusters
##############################
earlyearly.endoderm <- c('Gsc', 'Trh' , 'Otx2') # pg2
hepatic.progenitor <- c('Ttr', 'Hhex' , 'Tbx3')# pg2
ectoderm <- c('Gm2694' ,'Mir124-2hg')# pg2
blood <- c('Fam212a', 'Mgst3','Smim1','Blvrb') # Fig S1
cardiac <- c('Wnt2','Sh3bgr','Tnni1','Acta2')



sce
#class: SingleCellExperiment 
#dim: 20483 11039 


# markers.up <- findMarkers(sce, test="wilcox", # if wilcox test rather than t-test, get AUC rather than lfc
#                           groups=sce$celltype,  #lfc=1,
#                           direction="up") #, block=sce$sample)
# names(markers.up) # 36 cell types
# 
# save(markers.up, file='markers.up_wilcox_celltype_mesoderm.RData', compress=T)
# dim(markers.up[[1]]) # 22520    40
# ################# plot wilcox test ###################
# 
# library(pheatmap)
# load('markers.up_altas_celltype.RData')
# 
# pdf(file='heatmap_wilcox_upgene_altas_celltype_top5.pdf', height=15)
# for(i in names(markers.up))
# {
#   chosen <- i # "9"
#   
#   interesting.up <- markers.up[[chosen]]
#   best.set <- interesting.up[interesting.up$Top <= 5,]
#   AUCs <- getMarkerEffects(best.set, prefix="AUC")
#   pheatmap(AUCs, breaks=seq(-5, 5, length.out=101),
#            fontsize = 10, main= paste('Cluster',i))
# }
# dev.off()

markers.up <- findMarkers(sce, test="t", # if wilcox test rather than t-test, get AUC rather than lfc
                          groups=sce$subcelltype,  lfc=1,
                          direction="up") #, block=mnn.all$sample)
save(markers.up, file='markers.up.lfc1_subcelltype_16mesoderm.RData.RData')
lengths(markers.up)
# blood              cardiac.a              cardiac.b 
# 20483                  20483                  20483 
# cardiac.c          endothelial.a          endothelial.b 
# 20483                  20483                  20483 
# endothelial.c          endothelial.d extraembryonicMesoderm 
# 20483                  20483                  20483 
# mesodermProgenitors        mixedMesoderm.a        mixedMesoderm.b 
# 20483                  20483                  20483 
# pharyngealMesoderm   presomiticMesoderm.a   presomiticMesoderm.b 
# 20483                  20483                  20483 
# somiticMesoderm 
# 20483
library(openxlsx)
#library(xlsx)


pdf(file='heatmap_lfc1.FDR0.01.top10_subcelltype_16mesoderm.pdf', height=7)
for(i in names(markers.up))
{
  chosen <- i # "9"
  
  interesting.up <- markers.up[[chosen]]
  best.set <- subset(interesting.up, summary.logFC > 1 & FDR< 0.01 & Top<=10)
  logFCs <- getMarkerEffects(best.set)
  # rownames(logFCs) <- rowData(dimred.sce)[rownames(logFCs),]$SYMBOL
  pheatmap(logFCs, breaks=seq(-5, 5, length.out=101),
           fontsize = 10, main= paste('Cluster',i))
  
  write.xlsx(as.data.frame(best.set), row.names=TRUE,
             file= paste0('subcelltype_16mesoderm_lfc1.FDR0.01.top10/cluster',chosen,'.xlsx')) 
  #       SheetName = paste0('cluster',chosen)
  #     } else {write.xlsx(as.data.frame(best.set), file= "Wilox.up.markers.xlsx", 
  #                        SheetName = paste0('cluster',chosen), append=TRUE)
  # }
}
dev.off()




######################################
## output R object to run BioTIP
######################################

load(file='F:/projects/scRNA/results/IbarraSoria2018_MouseE8.25/normalized.RData')
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
lengths(samplesL)
# extraembryonicMesoderm          endothelial.a          endothelial.b          endothelial.c
#            1017                    239                    128                    458
# endothelial.d                  blood    mesodermProgenitors   presomiticMesoderm.b
#           47                   2079                   1696                    275
# presomiticMesoderm.a        somiticMesoderm        mixedMesoderm.a     pharyngealMesoderm
#         875                    739                   1328                    743
# mixedMesoderm.b              cardiac.a              cardiac.b              cardiac.c
#         715                    234                    184                    282

sce <- sce[,unlist(samplesL)]
sce
# dim: 20483 11039
sce$celltype <- as.factor(as.vector(sce$celltype))
sce$subcelltype <- as.factor(as.vector(sce$subcelltype))
sce$label <- as.factor(as.vector(sce$label))


## select global VG ---------------------------------
##################################

library(scran)
dec <- modelGeneVar(sce)

rownames(sce) <- uniquifyFeatureNames(rowData(sce)$GENEID, rowData(sce)$SYMBOL)

# Ordering by most interesting genes for inspection.
dec[order(dec$bio, decreasing=TRUE),] 

chosen <- getTopHVGs(dec, n=4000)
str(chosen)
length(chosen) # 4000
length(unique(chosen)) # 4000

sce <- sce[chosen, ]

save(sce, 
     file='F:/projects/scRNA/results/IbarraSoria2018_MouseE8.25/E8.25.2018_robustness/sce_16subtype.RData', 
     compress=TRUE)


#############################################
## output for files to run QuanTC
#############################################
{
  ## cell-cell similarity matrix for 131 cells #---------------------
  # refer to http://127.0.0.1:11637/library/SC3/doc/SC3.html
  # needs to customize ks accordingily per dataset, the larger range the longer running time.
  # In this case, ks=9:19 are tested, 
  # and for the related soft-thresholding clustering (QuanTC method), 
  # we had take average of the Consensus.Cluster-agreeable clustering results of k=4:10 to get cell-cell similarity matrix M
  library(SC3)
  sce.sc3 = sce
  rowData(sce.sc3)$feature_symbol <- rownames(sce.sc3)
  # remove features with duplicated names
  if(any(duplicated(rowData(sce.sc3)$feature_symbol))) sce.sc3 <- sce.sc3[-which(duplicated(rowData(sce.sc3)$feature_symbol)),]  # F
  range(assays(sce.sc3)$logcounts)
  # [1]   0.00000 11.34346
  
  ### to run SC3 successfully, transform sparsematrix to matrix  !!!!!!!!!!!!!
  logcounts(sce.sc3) <- as.matrix(logcounts(sce.sc3))
  counts(sce.sc3) <- as.matrix(counts(sce.sc3))
  
  # NOT repeat !!!! 
  set.seed(2020)
  sce.sc3 <- sc3(sce.sc3, ks = c(5,10, 14,16,18, 20), biology = FALSE) # svm_max = 5000 is default!!!
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
  # $ k_estimation   : num 31
  
  # to save space, transform back matrix to sparse matrix
  assayNames(sce.sc3)
  #[1]  "counts" "logcounts"
  save(sce.sc3, file='E8.25.2018_robustness/sce_SC3.RData', compress=TRUE) 
  assays(sce.sc3) <- assays(sce.sc3)[1]
  
  save(sce.sc3, file='E8.25.2018_robustness/sce_SC3.RData', compress=TRUE) 
  gc()
  # END DO NOT REPET !!!!!!!!!!!!!
  
  ##  writing input fiels for QuanTC -----------------------------------------
 # M_5 = (sce.sc3@metadata$sc3$consensus$`5`$consensus)
  M_10 = (sce.sc3@metadata$sc3$consensus$`10`$consensus)
  M_16 = (sce.sc3@metadata$sc3$consensus$`16`$consensus)
  M_20 = (sce.sc3@metadata$sc3$consensus$`20`$consensus)
  # take average of the Consensus.Cluster-agreeable clustering results to get cell-cell similarity matrix M
  M = (M_10+M_16+M_20)/3 
  
  QuanTC_input_dir ='F:/projects/QuanTC/QuanTC-modified/Input/E8.25.2018/'
  write.table(M, file=paste0(QuanTC_input_dir,'cell-cell.csv'), 
              row.names=FALSE, col.names = FALSE, sep=',')   
  
  logmat <- as.matrix(logcounts(sce))
  dim(logmat)
  # [1] 4000 11039
  logmat <- logmat[rownames(sce.sc3),]
  dim(logmat)
  # [1] 3999 1531
  write.table(logmat, file=paste0(QuanTC_input_dir, 'logmat.txt'), 
              sep='\t', row.names=FALSE, col.names=FALSE)
  
  true_label = as.vector(colData(sce)$subcelltype)
  write.table(true_label, file=paste0(QuanTC_input_dir, 'true_label.txt'), 
              sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  gene_name = rownames(logmat)
  length(gene_name)
  #[1] 3999
  gene_name[1:4]
  #[1]"Hbb-bh1" "Hba-x"   "Hba-a1"  "Phlda2" 
  write.table(gene_name, file=paste0(QuanTC_input_dir, 'gene_name.txt'), 
              sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE)
  
}


