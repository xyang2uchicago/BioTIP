
setwd("E:/Git_Holly/scRNAseq_examples/result/simulated_EMT")

library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(BioTIP)


####################################################################
## generate the Ic scores for 100 random genes
####################################################################
load('../../data/simulated_EMT/sce.RData')
sce
# class: SingleCellExperiment 
# dim: 18 5363 
# metadata(1): study
# assays(2): counts logcounts
# rownames(18): snailt SNAIL ... Ncad Ovol2
# rowData names(1): gene_name
# colnames(5363): 1 2 ... 5362 5363
# colData names(8): true_label cell_type ... C_Leiden_0.8 C_Leiden_1.2
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):
  
table(sce$cell_type)
# E   I1   I2    M   TC 
# 577 1849  801  985 1151

logmat <- as.matrix(logcounts(sce))
samplesL <- split(colnames(logmat), f = factor(sce$cell_type))
lengths(samplesL)
# E   I1   I2    M   TC 
# 577 1849  801  985 1151

n.sim = 100
n.gene = 10
Ic <- matrix(nrow=n.sim, ncol=length(samplesL))

#set.seed(2022)
for(i in 1:n.sim){
  random.gene <- sample(rownames(logmat), n.gene)
  Ic[i,] <- getIc(logmat, samplesL, random.gene, fun="cor")
}
colnames(Ic) <- names(samplesL)
save(Ic, file=paste0('existing.Ic_',n.gene,'gene.RData') ) 

# reorder according to pesudo orders
Ic <- Ic[,c('E',   'I1' ,   'TC',  'I2',    'M'  )]
         
average.Ic <- apply(Ic,2,mean)

pdf(file='average.Ic_10randomgene.pdf', height=3)
boxplot(Ic,xaxt='n', 
        xlab='cell cluster', main='EMT',
        las=1)
lines(1:length(average.Ic), average.Ic, type='b', col='red')
axis(1, at=1:length(average.Ic),labels= names(average.Ic), las=2)
dev.off()

wilcox.test(Ic[,'I1'], Ic[,'TC']) # p-value < 2.2e-16
wilcox.test(Ic[,'TC'], Ic[,'I2'])  # p-value < 2.2e-16



 

