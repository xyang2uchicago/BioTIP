
setwd("E:/Git_Holly/scRNAseq_examples/result/gastrulationE8.25_Pijuan-Sala2019")

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
load('../../data/gastrulationE8.25_Pijuan-Sala2019/sce_E8.25_HEP.RData')
sce
# class: SingleCellExperiment 
# dim: 3073 1362 
# metadata(0):
# assays(2): counts logcounts
# rownames(3073): Phlda2 Myl7 ... Hmcn1 Tfdp2
# rowData names(2): ENSEMBL SYMBOL
# colnames(1362): cell_63240 cell_63242 ... cell_95715 cell_95717
# colData names(31): cell barcode ... C_Leiden_0.8 C_Leiden_1.2
# reducedDimNames(5): pca.corrected umap TSNE UMAP force
# altExpNames(0):

logmat <- as.matrix(logcounts(sce))
samplesL <- split(rownames(colData(sce)), f =sce$label)
all(table(sce$label) == lengths(samplesL))  #TRUE

n.sim = 100
Ic <- matrix(nrow=n.sim, ncol=length(samplesL))

set.seed(2022)
for(i in 1:n.sim){
  random.gene <- sample(rownames(logmat), n.sim)
  Ic[i,] <- getIc(logmat, samplesL, random.gene, fun="cor")
}
colnames(Ic) <- paste0('C', names(samplesL))
save(Ic, file='existing.Ic_100gene.RData')  


average.Ic <- apply(Ic,2,mean)
pdf(file='average.Ic_100randomgene.pdf', height=3)
boxplot(Ic,xaxt='n', 
        xlab='cell cluster', main='E8.25 2019',
        las=1)
lines(1:length(average.Ic), average.Ic, type='b', col='red')
axis(1, at=1:length(average.Ic),labels=names(average.Ic), las=2)
dev.off()

wilcox.test(Ic[,'C6'], Ic[,'C15']) # p<2e-16
wilcox.test(Ic[,'C15'], Ic[,'C7'])  #  p<2e-16

wilcox.test(Ic[,'C3'], Ic[,'C13']) # p-value <2e-16
wilcox.test(Ic[,'C13'], Ic[,'C10']) # p-value = 0.02975






