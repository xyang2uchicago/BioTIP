
setwd('F:/projects/BioTIP/result/QuanTC_simulation/simulatedEMT_robustness/')
load(file='../sce.RData')

################################################
## replot simulated DNB scores with extended y-axis 
###########################################
library(BioTIP)
load('C_cell_type/BioTIP.res.RData')
load('C_cell_type/BioTIP_top1FDR0.05_SimuMCI_100_1_fdr0.05_minsize6.RData')
reorder <- c('E','I1','TC','I2','M')

pdf(file='C_cell_type/BioTIP_top1FDR0.05_barplot_MCI_Sim_RandomGene.pdf',height=3, width=9)
par(mfrow=c(1,3))
for(i in 1:3){
  plot_MCI_Simulation(res$CTS.score[i], simuMCI[[i]][reorder,], las=2, ylim=c(0, 0.5), na.rm=TRUE)
}
dev.off()


################################################
## replot simulated IC.shrink scores with extended y-axis 
###########################################
library(BioTIP)
load('C_cell_type/BioTIP.res.RData')
load('C_cell_type/BioTIP_top1FDR0.05_IC_sim.PermutateBoth.RData')
reorder <- c('E','I1','TC','I2','M')

nx <- length(res$CTS.candidate)
pdf(file='C_cell_type/BioTIP_top1FDR0.05_IC_Delta_SimresultBoth.pdf',height=5, width=7)

par(mfrow=c(2,nx))
for(i in 1:nx){
  x = length(res$CTS.candidate[[i]])
  plot_Ic_Simulation(BioTIP_scores[[i]][reorder], SimResults_b[[i]][reorder,], 
                     ylim=c(0, 1),
                     las = 2, ylab="Ic.shrink",
                     main= paste("Cluster",names(BioTIP_scores)[i],"\n",x," genes"),
                     fun="boxplot",  #fun="matplot", 
                     which2point=names(BioTIP_scores)[i])
 }
for(i in 1:nx){
  x = length(res$CTS.candidate[[i]])
  plot_SS_Simulation(BioTIP_scores[[i]][reorder], SimResults_b[[i]][reorder,], 
                     main = paste("Delta Ic.shrink",x,"genes"), 
                     ylab=NULL,
                     xlim=range(0,0.2))
}

dev.off()

###############################
## PLOT gene expression changes
################################
table(sce$C_SNNGraph.k100, sce$cell_type)
#      E   I1   I2    M   TC
# 1    0 1783    0    0   15
# 2    0   16  755   22   92
# 3    0   50   38    5  618
# 4  577    0    0    0  176
# 5    0    0    8  958  250
table(sce$C_SNNGraph.k200, sce$cell_type)
#      E   I1   I2    M   TC
# 1    0 1830   31    0  512
# 2    0   19  762   25  211
# 3  577    0    0    0  175
# 4    0    0    8  960  253

samplesL <- split(rownames(colData(sce)),f = colData(sce)$cell_type) 
lengths(samplesL)
#  E   I1   I2    M   TC 
#577 1849  801  985 1151  

logmat <- as.matrix(logcounts(sce))
tmp <- data.frame(t(logmat[, unlist(samplesL)]))
tmp$cell_type= factor(colData(sce[, unlist(samplesL)])$cell_type,
                      levels=c('E','I1','TC','I2','M'))
tmp$C_SNNGraph.k100= factor(colData(sce[, unlist(samplesL)])$C_SNNGraph.k100,
                            levels = c('4', '1',  '3',   '2', '5'))
# sce$correct_cluster = rep('F', ncol(sce))  
# sce$correct_cluster[which(as.vector(sce$cell_type)=='E' & sce$C_SNNGraph.k200=='4')] <- 'T'
# sce$correct_cluster[which(as.vector(sce$cell_type)=='I1' & sce$C_SNNGraph.k200=='1')] <- 'T'
# sce$correct_cluster[which(as.vector(sce$cell_type)=='I2' & sce$C_SNNGraph.k200=='2')] <- 'T'
# sce$correct_cluster[which(as.vector(sce$cell_type)=='M' & sce$C_SNNGraph.k200=='5')] <- 'T'
# sce$correct_cluster[which(as.vector(sce$cell_type)=='TC' & sce$C_SNNGraph.k200=='3')] <- 'T'
tmp$C_SNNGraph.k200= factor(colData(sce[, unlist(samplesL)])$C_SNNGraph.k200,
                           levels = c('3', '1',  '2',   '4'))
sce$correct_cluster = rep('F', ncol(sce))  
sce$correct_cluster[which(as.vector(sce$cell_type)=='E' & sce$C_SNNGraph.k200=='3')] <- 'T'
sce$correct_cluster[which(as.vector(sce$cell_type)=='I1' & sce$C_SNNGraph.k200=='1')] <- 'T'
sce$correct_cluster[which(as.vector(sce$cell_type)=='I2' & sce$C_SNNGraph.k200=='2')] <- 'T'
sce$correct_cluster[which(as.vector(sce$cell_type)=='M' & sce$C_SNNGraph.k200=='4')] <- 'T'
tmp$correct_cluster <- factor(sce[, unlist(samplesL)]$correct_cluster,
                              levels=c('T','F'))

levels(sce$cell_type) <- c("E" , "I1", "TC", "I2", "M"  )
mycolor = c( "dodgerblue1", "orange", "green","olivedrab2","purple")
names(mycolor) <- levels(sce$cell_type)

pdf(file='boxplot_10genes.pdf', width=10, height=4) 
gridExtra::grid.arrange(
  ggplot(tmp, aes(x=cell_type, y=miR200t, color=cell_type)) + geom_boxplot(show.legend = FALSE) +
    scale_fill_manual(values=mycolor),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=tgfR, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=tgft, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=ZEB, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=zebt, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=ZR1, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=ZR2, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=ZR3, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=ZR4, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=cell_type, y=ZR5, color=cell_type)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ncol=5)
gridExtra::grid.arrange(
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=miR200t, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=tgfR, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=tgft, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=ZEB, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=zebt, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=ZR1, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=ZR2, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=ZR3, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=ZR4, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ggplot(tmp, aes(x=C_SNNGraph.k200, y=ZR5, color=correct_cluster)) + geom_boxplot(show.legend = FALSE),  #geom_point(),
  ncol=5)
dev.off()

pdf('umap.pdf')
plotReducedDim(sce, dimred="UMAP", colour_by='cell_type', 
               text_by='cell_type',
               point_size=0.2, text_siz=4) 
plotReducedDim(sce, dimred="UMAP", colour_by='cell_type', 
               text_by='C_SNNGraph.k200',
               point_size=0.2, text_siz=4) 
plotReducedDim(sce, dimred="UMAP", colour_by='cell_type', 
               text_by='C_SNNGraph.k100',
               point_size=0.2, text_siz=4) 
dev.off()
