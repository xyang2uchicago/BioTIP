setwd('F:/projects/scRNA/results/GSE87038_gastrulation/uncorrected/E8.25_mesoderm_robustness')

# the following funciton has been updated on github, but before updating the package to Bioconductor, you need to call it locally to test
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/ForJennifer/optimize.sd_selection.R') 
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/BioTIP.wrap.R')
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/stability_score.R') 

# normalized data were downloaded and reanalyzed
# refer to GSE87038_E8.25_normalize.R
library(BioTIP)
library(dplyr)
library(SingleCellExperiment)

library(ggplot2)
library(gridExtra)  # required by grid.arrange()
library(ggpubr)  # required by p-values on bar

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

table(sce$subcelltype)
# amnion.a               amnion.b                  blood              cardiac.a
# 237                    533                   2079                    234
# cardiac.b              cardiac.c          endothelial.a          endothelial.b
# 184                    282                    239                    128
# endothelial.c          endothelial.d extraembryonicEctoderm extraembryonicEndoderm
# 458                     47                    367                   1876
# extraembryonicMesoderm              forebrain                foregut    mesodermProgenitors
# 1017                    764                    185                   1696
# midHindbrain           midHindgut.a           midHindgut.b        mixedMesoderm.a
# 1242                    158                     95                   1328
# mixedMesoderm.b            neuralCrest             neuralTube            notochord.a
# 715                    332                   1128                     41
# notochord.b            notochord.c     pharyngealMesoderm             placodes.a
# 17                     27                    743                    805
# placodes.b             placodes.c   presomiticMesoderm.a   presomiticMesoderm.b
# 439                    101                    875                    275
# somiticMesoderm
# 739

X <- logcounts(sce)
dim(X) # 30483  19386
y <- apply(X,1,sum, na.rm=T)
table(y>0)
# FALSE  TRUE
# 133   20350
X <- X[y>0,]
dim(X)  #[1] 20350 19386

logmat <- as.matrix(X)
dim(logmat) # 20350 19386



load('original_run/CTS.RData')
lengths(CTS)
# 15  6 16 13 
# 67 90 79 60
# formating gene alias with official symbols
CTS.Lib.Symbol <- CTS
CTS.Lib.Symbol[[1]][which(CTS.Lib.Symbol[[1]] == 'Sept1')] = 'Septin1'
CTS.Lib.Symbol[[1]][which(CTS.Lib.Symbol[[1]] == 'Grasp')] = 'Tamalin'
CTS.Lib.Symbol[[1]][which(CTS.Lib.Symbol[[1]] == 'Fam57b')] = 'Tlcd3b'
CTS.Lib.Symbol[[3]][which(CTS.Lib.Symbol[[3]] == 'Sept4')] = 'Septin4'
CTS.Lib.Symbol[[3]][which(CTS.Lib.Symbol[[3]] == 'Fam69a')] = 'Dipk1a'
CTS.Lib.Symbol[[3]][which(CTS.Lib.Symbol[[3]] == 'Sdpr')] = 'Cavin2'
CTS.Lib.Symbol[[3]][which(CTS.Lib.Symbol[[3]] == 'Wisp1')] = 'Ccn4'
CTS.Lib.Symbol[[4]][which(CTS.Lib.Symbol[[4]] == 'Cxx1b')] = 'Rtl8b'

setdiff(CTS.Lib.Symbol[[1]], rownames(logmat)) # 0    
setdiff(CTS.Lib.Symbol[[2]], rownames(logmat)) # 0
setdiff(CTS.Lib.Symbol[[3]], rownames(logmat))  # 0
setdiff(CTS.Lib.Symbol[[4]], rownames(logmat)) # 0

#par(mfrow=c(2,2)) 
for(i in 1:4)
{
  # testing for those measured CTS gene members
  CTS <-intersect(CTS.Lib.Symbol[[i]], rownames(logmat))
  IC.shrink <- getIc(logmat[,unlist(samplesL)], samplesL, CTS, fun="BioTIP",
                         shrink=TRUE, PCC_sample.target='none'
  )
  state.in.order = names(samplesL)
  n <- length(IC.shrink)
  par(mar=c( 12,5,5,5))
  plot(IC.shrink, type='b',xaxt ="n", main=paste0('CTS: E8.25 C', names(CTS.Lib.Symbol)[i]))
  axis(side=1, at= 1:n, labels=state.in.order,  cex.axis=1, las=2)
}


load('../E8.25_mesoderm/BioTIP_top0.1FDR0.2/IC_sim_S13_CTS60_PermutateGene.in.IbarraSoria2018_PCCsNoShrink.RData')
