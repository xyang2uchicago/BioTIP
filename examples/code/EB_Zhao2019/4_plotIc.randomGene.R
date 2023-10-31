
#setwd("F:/projects/scRNA/results/AJ/GSE130146_xy/Results_1k/LibrarySize/GSE130146_robustness/") 
setwd("E:/Git_Holly/scRNAseq_examples/result/EB_Zhao2019")

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
load('../../data/EB_Zhao2019/sce.GSE130146_noenderdormPgerm.RData')
sce
# class: SingleCellExperiment 
# dim: 4000 1531 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(4000): Mesp1 Dppa5a ... 6330408A02Rik Tmem231
# rowData names(2): ID Symbol
# colnames(1531): AAACCTGAGTTTCCTT-1 AAACCTGCATGAAGTA-1 ...
# TTTGTCATCGGCTTGG-1 TTTGTCATCTGCTTGC-1
# colData names(14): Sample Barcode ... C_Leiden_1.2 C_Soft
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

table(sce$C_SNNGraph_k5)
#   1  10  11  12  14  17   3   4   5   6   7   8   9 
# 102  38 208  34  39  11 289  98  86  75 151 160 240
logmat <- as.matrix(logcounts(sce))
samplesL <- split(rownames(colData(sce)), f = sce$C_SNNGraph_k5)
lengths(samplesL)

n.sim = 100
n.gene = 50
Ic <- matrix(nrow=n.sim, ncol=length(samplesL))
#  1  10  11  12  14  17   3   4   5   6   7   8   9 
# 102  38 208  34  39  11 289  98  86  75 151 160 240 

set.seed(2022)
for(i in 1:n.sim){
  random.gene <- sample(rownames(logmat), n.gene)
  Ic[i,] <- getIc(logmat, samplesL, random.gene, fun="cor")
}
colnames(Ic) <- paste0('C',names(samplesL))
save(Ic, file=paste0('existing.Ic_',n.gene,'gene.RData') ) 

# reorder according to pesudo orders
Ic <- Ic[,c('C12', 'C6','C7','C10', 'C9', 'C4','C1','C3', 'C11', 'C8', 'C5', 'C14', 'C17' )]
average.Ic <- apply(Ic,2,mean)

pdf(file='average.Ic_50randomgene.pdf', height=3)
boxplot(Ic,xaxt='n', 
        xlab='cell cluster', main='EB',
        las=1)
lines(1:length(average.Ic), average.Ic, type='b', col='red')
axis(1, at=1:length(average.Ic),labels= names(average.Ic), las=2)
dev.off()


wilcox.test(Ic[,'C9'], Ic[,'C4']) # p-value < 2.2e-16
wilcox.test(Ic[,'C4'], Ic[,'C1']) # p-value =  0.004505


