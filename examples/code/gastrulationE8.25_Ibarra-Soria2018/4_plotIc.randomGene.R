
setwd("E:/Git_Holly/scRNAseq_examples/result/gastrulationE8.25_Ibarra-Soria2018")

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
load('../../data/gastrulationE8.25_Ibarra-Soria2018/sce_16subtype.RData')
sce
# class: SingleCellExperiment 
# dim: 4000 11039 
# metadata(0):
#   assays(3): counts normcounts logcounts
# rownames(4000): Hbb-bh1 Hba-x ... 2900060B14Rik Phf8
# rowData names(5): GENEID SYMBOL SEQNAME HVGs.10 HVGs.20
# colnames(11039): extraembryonicMesoderm_1 extraembryonicMesoderm_2 ...
# cardiac.c_281 cardiac.c_282
# colData names(3): label celltype subcelltype
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

logmat <- as.matrix(logcounts(sce))
samplesL <- split(rownames(colData(sce)), f =sce$subcelltype)
all(lengths(samplesL)==table(sce$subcelltype)) #TRUE

n.sim = 100
n.gene = 50
Ic <- matrix(nrow=n.sim, ncol=length(samplesL))

set.seed(2022)
for(i in 1:n.sim){
  random.gene <- sample(rownames(logmat), n.gene)
  Ic[i,] <- getIc(logmat, samplesL, random.gene, fun="cor")
}
colnames(Ic) <- names(samplesL)
save(Ic, file='existing.Ic_50gene.RData')  


average.Ic <- apply(Ic,2,mean)
max(average.Ic)

pdf(file='average.Ic_50randomgene.pdf', height=3)
boxplot(Ic,xaxt='n', 
        xlab='cell subtype', main='EB',
        las=1)
lines(1:length(average.Ic), average.Ic, type='b', col='red')
axis(1, at=1:length(average.Ic),labels=names(average.Ic), las=2)
dev.off()

wilcox.test(Ic[, "endothelial.b"], Ic[, "endothelial.a"]) #p-value < 2.2e-16
wilcox.test(Ic[, "endothelial.b"], Ic[, "endothelial.c"]) #p-value < 2.2e-16

wilcox.test(Ic[, "blood"], Ic[, "endothelial.d"]) #p-value < 2.2e-16


