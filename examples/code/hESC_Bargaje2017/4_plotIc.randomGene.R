
setwd("E:/Git_Holly/scRNAseq_examples/result/hESC_Bargaje2017")

# PubMed ID	29311656
# GRCm38.p4
# annotation from the Ensembl database, version 84.
#
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(BioTIP)


####################################################################
## generate the Ic scores for 100 random genes
####################################################################
load('../../data/hESC_Bargaje2017/sce_hESC.RData')
sce
# class: SingleCellExperiment 
# dim: 96 929 
# metadata(0):
#   assays(1): logcounts
# rownames(96): ACVR1B ACVR2A ... WNT5A WNT5B
# rowData names(0):
#   colnames(929): Cell1 Cell109 ... Cell1421 Cell1429
# colData names(17): CollectionTime Consensus.Cluster ... C_Leiden_0.8
# C_Leiden_1.2
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

table(sce$Consensus.Cluster) 
# 1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
# 77 107   0   0   0   0   0   0   0   0 161  86  86  80  15 173  95  49 
logmat <- as.matrix(logcounts(sce))
samplesL <- split(rownames(colData(sce)), f =sce$Consensus.Cluster)
lengths(samplesL)
# 1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
# 77 107   0   0   0   0   0   0   0   0 161  86  86  80  15 173  95  49 
samplesL = samplesL[which(lengths(samplesL)>20)]

n.sim = 100
n.gene = 20
Ic <- matrix(nrow=n.sim, ncol=length(samplesL))

set.seed(2022)
for(i in 1:n.sim){
  random.gene <- sample(rownames(logmat), n.gene)
  Ic[i,] <- getIc(logmat, samplesL, random.gene, fun="cor")
}
colnames(Ic) <- paste0('C',names(samplesL))
save(Ic, file='existing.Ic_20gene.RData')  

# reorder according to pesudo orders
Ic <- Ic[,c( "C2", "C1" , "C3",  "C4",  "C5", "C9", "C7" , "C8","C10"  )]

average.Ic <- apply(Ic,2,mean)

pdf(file='average.Ic_20randomgene.pdf', height=3)
boxplot(Ic,xaxt='n', 
        xlab='cell cluster', main='hESC_Bargaje2017',
        las=1)
lines(1:length(average.Ic), average.Ic, type='b', col='red')
axis(1, at=1:length(average.Ic),labels= names(average.Ic), las=2)
dev.off()


wilcox.test(Ic[,'C9'], Ic[,'C5']) # p= 5.286e-05
wilcox.test(Ic[,'C9'], Ic[,'C7'])  #  p<2e-16

wilcox.test(Ic[,'C8'], Ic[,'C10']) # p-value = 0.1686

 



