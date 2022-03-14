
setwd("E:/Git_Holly/scRNAseq_examples/result/lung_Treutlein2014")

library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(BioTIP)


####################################################################
## generate the Ic scores for 100 random genes
####################################################################
load('../../data/lung_Treutlein2014/AT2.sce.RData')
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

logmat <- as.matrix(assay(sce, 'logFPKM'))
samplesL <- split(rownames(colData(sce)), f = sce$C_Leiden_0.4)
table(sce$C_Leiden_0.4)
# 1  2  3  4  5 
# 46 32 24 17 12 
lengths(samplesL)
# 1  2  3  4  5 
# 46 32 24 17 12

n.sim = 100
n.gene = 50
Ic <- matrix(nrow=n.sim, ncol=length(samplesL))

#set.seed(2022)
for(i in 1:n.sim){
  random.gene <- sample(rownames(logmat), n.gene)
  Ic[i,] <- getIc(logmat, samplesL, random.gene, fun="cor")
}
colnames(Ic) <- paste0('C',names(samplesL))
save(Ic, file=paste0('existing.Ic_',n.gene,'gene.RData') ) 

# reorder according to pesudo orders
Ic <- Ic[,c('C3','C4','C2','C5','C1' )]
average.Ic <- apply(Ic,2,mean)

pdf(file='average.Ic_50randomgene.pdf', height=3)
boxplot(Ic,xaxt='n', 
        xlab='cell cluster', main='lung',
        las=1)
lines(1:length(average.Ic), average.Ic, type='b', col='red')
axis(1, at=1:length(average.Ic),labels= names(average.Ic), las=2)
dev.off()


wilcox.test(Ic[,'C2'], Ic[,'C1']) # p-value 0.0002346
wilcox.test(Ic[,'C2'], Ic[,'C4'])  # p-value =0.1842


