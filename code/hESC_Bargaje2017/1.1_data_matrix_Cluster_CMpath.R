# R4.0.2
 setwd('F:/projects/BioTIP/result/Bargaje2017_EOMES/Cluster')

 library(RColorBrewer)
 library(SingleCellExperiment)
 library(org.Mm.eg.db)
 library(scater)
 library(scran)
 library(pheatmap)
 library(limma)
 library(edgeR)
 library(dynamicTreeCut)
 library(cluster)
 library('cluster')
 library(bluster)
 
 
###################################################
## 1) read the published data matrix of 96 genes ##
###################################################

###########################################
## 2) reduce dimention and visualization ##
###########################################

################################################################################
## 3) defining cell identity by the expression levels of known lineage-marker ##
################################################################################
# here we used the published consensus clusters
 
###################################################
## 4) trajectory analysis ##
###################################################
 load('../sce.RData')  
{
  sec
 # class: SingleCellExperiment 
 # dim: 96 1896 
 
 table(sce$Consensus.Cluster)
 #  1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
 # 77 107 130 299 157  64  70  81  67  77 161  86  86  80  37 173  95  49
 
 tmp <- subset(colData(sce), CollectionTime=='Day3' )
 int <- which(colData(sce)$CollectionTime %in% c('Day0', 'Day1', 'Day1.5', 'Day2', 'Day2.5') |
                (colData(sce)$CollectionTime =='Day3' & colData(sce)$Consensus.Cluster ==10))
 length(int)  #[1] 929
 sce <- sce[,int]
 sce$Consensus.Cluster <- factor(as.vector(sce$Consensus.Cluster))
 table(sce$Consensus.Cluster)
 # 1  10   2   3   4   5   6   7   8   9
 # 77 107 161  86  86  80  15 173  95  49

 
 library(scater)
 # pseudo counts
 #counts(sce) = 2^(logcounts(sce))
 by.cluster <- aggregateAcrossCells(sce, ids=sce$Consensus.Cluster, use.assay.type='logcounts')
 centroids <- reducedDim(by.cluster, "PCA")

 dmat <- dist(centroids)
 dmat <- as.matrix(dmat)
 g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
 mst <- igraph::minimum.spanning.tree(g)
 
 pdf(file="../hESC_Bargaje2017_robustness/trajectory_929cells.pdf")

 #set.seed(1000)
 plot(mst)
 
 pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
 coords <- reducedDim(by.cluster, "TSNE")
 group <- rep(seq_len(nrow(pairs)), 2)

 stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)
 
 plotTSNE(sce, colour_by="Consensus.Cluster",
          text_by="Consensus.Cluster", text_size = 8)
 
 plotTSNE(sce, colour_by="Consensus.Cluster", 
          text_by="Consensus.Cluster", text_size = 8) + 
   geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
 
 
 dev.off()
} 

 
 
 
############################
## 4) BioTIP application  ##
############################
library(BioTIP)
packageVersion('BioTIP') # '1.5.0'

# Ic, we generated time point-specific Log2Ex matrices taht were downloaded from teh publication 
# For each of them (days 0, 1, 1.5, 2, 2.5, and 3/M-only mesoderm-specific cells)
tmp <- subset(colData(sce), CollectionTime=='Day3' )
table(tmp$Consensus.Cluster)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18 
# 0   0   0   0   0  22   0   0   0 106 130   0   0   0   0   0   0   0 

int <- which(colData(sce)$CollectionTime %in% c('Day0', 'Day1', 'Day1.5', 'Day2', 'Day2.5') |
               (colData(sce)$CollectionTime =='Day3' & colData(sce)$Consensus.Cluster ==10))
length(int)  #[1] 929

table(colData(sce)[int,]$CollectionTime)
#  Day0   Day1 Day1.5   Day2 Day2.5   Day3   Day4   Day5 
#  231    166     93    211    122    106      0      0 


sce <- sce[,int]
sce
# class: SingleCellExperiment 
# dim: 96 929 
# metadata(0):
#   assays(1): logcounts
# rownames(96): ACVR1B ACVR2A ... WNT5A WNT5B
# rowData names(0):
#   colnames(929): Cell1 Cell109 ... Cell1421 Cell1429
# colData names(3): CollectionTime Consensus.Cluster Phenotype.FACS
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):
  
# Note that BioTIP further ignores 15 Day2-2.5 cells that are Endodern cluster 6 because its population size is <20!!
table(sce$CollectionTime,sce$Consensus.Cluster)
#          1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9
# Day0    74   0   0   0   0   0   0   0   0   0 154   3   0   0   0   0   0   0
# Day1     2   0   0   0   0   0   0   0   0   0   3  75  85   1   0   0   0   0
# Day1.5   1   0   0   0   0   0   0   0   0   0   4   8   1  79   0   0   0   0
# Day2     0   0   0   0   0   0   0   0   0   0   0   0   0   0   7 131  42  31
# Day2.5   0   1   0   0   0   0   0   0   0   0   0   0   0   0   8  42  53  18
# Day3     0 106   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# Day4     0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# Day5     0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0

#now convert year to a character and type to a factor
df <- table(sce$CollectionTime,sce$Consensus.Cluster)
df2 <- NULL
for(i in 1:ncol(df)) {
  tmp <- data.frame(Freq=df[,i], 
                    Cluster= paste0('C',rep(colnames(df)[i], nrow(df))),
                    Time=rownames(df)
  )
  df2 <- rbind(df2, tmp)
  }

df2 <- subset(df2, Freq>0)
head(df2)
#         Freq Cluster   Time
# Day0      74       C1   Day0
# Day1       2       C1   Day1
# Day1.5     1       C1 Day1.5
# Day2.51    1      C10 Day2.5
# Day31    106      C10   Day3
df2$Time <- factor(df2$Time, levels=c("Day0",   "Day1",   "Day1.5","Day2", "Day2.5", "Day3"))

# Calculate the cumulative sum of len for each dose
library(plyr)
df_cumsum <- ddply(df2, "Time",
                   transform, label_ypos=cumsum(Freq))
head(df_cumsum)
#    Freq Cluster Time label_ypos
# 1   74      C1 Day0         74
# 2  154      C2 Day0        228
# 3    3      C3 Day0        231
# 4    2      C1 Day1          2

ggplot(data=df_cumsum, aes(x=Time, y=Freq, fill=Cluster)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=Cluster), vjust=1.6, 
            color="white", size=3.5)
  
dev.copy2pdf(file='bar_Time_vs_CLuster.pdf')

int <- which(colData(sce)$Consensus.Cluster !=6)  # because it has only 15 cells
length(int)  #[1] 914

table(colData(sce)[int,]$CollectionTime)
# Day0   Day1 Day1.5   Day2 Day2.5   Day3   Day4   Day5 
# 231    166     93    204    114    106      0      0 

sce <- sce[,int]

# 1) module construction

# Pre-selection Transcript
cut.preselect = 0.8  # out of 96 genes
cut.fdr = 0.2   
cut.minsize = 10    

samplesL <- split(rownames(colData(sce)),f = colData(sce)$Consensus.Cluster) #!!!!!!!!!!!!
lengths(samplesL)
#  1  10  11  12  13  14  15  16  17  18   2   3   4   5   6   7   8   9 
# 77 107   0   0   0   0   0   0   0   0 161  86  86  80   0 173  95  49 
samplesL <- samplesL[-which(lengths(samplesL)<10)]
lengths(samplesL)
#  1  10   2   3   4   5   7   8   9 
# 77 107 161  86  86  80 173  95  49

logmat <- as.matrix(logcounts(sce))
dim(logmat) # [1] 96 914

#this optimizes sd selection
set.seed(2020)
testres <- optimize.sd_selection(logmat, samplesL, B=100, cutoff=cut.preselect,
                                 times=.75, percent=0.8)
class(testres) #[1] "list"
class(testres[[1]])  #[1] "matrix" "array" 
lapply(testres, dim)

# Network Partition
igraphL <- getNetwork(testres, fdr = cut.fdr)
# 1:63 nodes
# 10:76 nodes
# 2:72 nodes
# 3:74 nodes
# 4:70 nodes
# 5:72 nodes
# 7:73 nodes
# 8:75 nodes
# 9:69 nodes


cluster <- getCluster_methods(igraphL)

##### plot network #################
####################################
names(igraphL)
#[1] "1"  "10" "2"  "3"  "4"  "5"  "7"  "8"  "9" 

library('igraph')
pdf(file='Network.view_C9_C10.pdf', width=10)
par(mfrow=c(1,2))
tmp = igraphL[["9"]]
E(tmp)$width <- E(tmp)$weight*3
V(tmp)$community= cluster[["9"]]$membership
mark.groups = table(cluster[["9"]]$membership)
#mark.groups = mark.groups[which(mark.groups>=10)]
colrs = rainbow(length(mark.groups), alpha = 0.3)
V(tmp)$label <- NA
plot(tmp, vertex.color=colrs[V(tmp)$community], vertex.size = 5,
     mark.groups=cluster[["9"]])
legend(1,1, paste0(names(mark.groups),sep=":",mark.groups), text.col=rainbow(length(mark.groups)))

tmp = igraphL[["10"]]
E(tmp)$width <- E(tmp)$weight*3
V(tmp)$community= cluster[["10"]]$membership
mark.groups = table(cluster[["10"]]$membership)
colrs = rainbow(length(mark.groups), alpha = 0.3)
V(tmp)$label <- NA
plot(tmp, vertex.color=colrs[V(tmp)$community], vertex.size = 5,
     mark.groups=cluster[["10"]])
legend(1,1, paste0(names(mark.groups),sep=":",mark.groups), text.col=rainbow(length(mark.groups)))

dev.off()

# 2.1) Identifying CTS, new MCI score #########
membersL <- getMCI(cluster,testres, adjust.size = FALSE, fun='BioTIP')
names(membersL)
#[1] "members" "MCI"     "sd"      "PCC"     "PCCo"  
save(membersL, file="membersL.RData")

pdf(file=paste0("MCIBar_", cut.preselect,
                "_fdr",cut.fdr,"_minsize",cut.minsize,".pdf"),
    width=10, height=5)
plotBar_MCI(membersL, ylim=c(0,10), minsize = 30)
plotBar_MCI(membersL, ylim=c(0,10), minsize = 20)
plotBar_MCI(membersL, ylim=c(0,10), minsize = cut.minsize)
plotBar_MCI(membersL, ylim=c(0,10), minsize = 5)

dev.off()
#the numbers in the parenthesis: they are total number of clusters, no control of the cell (sample) size per  cluster


# get the statistics using the MCI system
n.state.candidate = 7
topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], 
                   min=cut.minsize, n=n.state.candidate)
names(topMCI)
#[1] "9"  "10" "7"  "1"  "4"  "8"  "5" 

# get the state ID and MCI statistics for the two leading MCI scores per state
maxMCIms <- getMaxMCImember(membersL[["members"]],
                            membersL[["MCI"]],min =cut.minsize, n=2)
names(maxMCIms)
#[1] "idx"             "members"         "2topest.members"  

maxMCI = getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
unlist(maxMCI)
#       1       10        2        3        4        5        7        8        9 
# 4.151014 5.276104 3.423058 3.002783 3.751412 3.583947 4.649539 3.640540 5.779724 

names(maxMCIms[["members"]][names(topMCI)])
CTS = getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
# Length: 18
# Length: 19
# Length: 10
# Length: 19
# Length: 23
# Length: 16
# Length: 27

## tmp calculates the number of bars within each named state
(tmp = unlist(lapply(maxMCIms[['idx']][names(topMCI)], length)))
# 9 10  7  1  4  8  5 
# 2  2  2  2  1  2  2

## here returns all the groups with exactly 2 bars
(whoistop2nd = names(tmp[tmp==2]))
#[1]  "9"  "10" "7"  "1"  "8"  "5" 

## add the gene members of the 2nd toppest biomodue in the states with exactly 2 bars
if(length(whoistop2nd)>0)  CTS = append(CTS, maxMCIms[["2topest.members"]][whoistop2nd])
names(CTS)
#[1][1] "9"  "10" "7"  "1"  "4"  "8"  "5"  "9"  "10" "7"  "1"  "8"  "5" 

# ## add the gene members of the 2nd toppest biomodue in the states with exactly 3 bars
# if(length(whoistop3rd)>0)  CTS = append(CTS, maxMCIms[["2topest.members"]][whoistop3rd])  
# names(CTS)
# #[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2"   "Day1"   "Day2.5" "Day1.5"

# ## add the gene members of the 3rd toppest biomodue in the states with exactly 3 bars
# if(length(whoistop3rd)>0)  CTS = append(CTS, maxMCIms[["3topest.members"]][whoistop3rd])  
# names(CTS)
# #[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2"   "Day1"   "Day2.5" "Day1.5" "Day2.5" "Day1.5"


'EOMES' %in% unlist(CTS)
#[1] FALSE
'SOX17' %in% unlist(CTS)
#[1] TRUE
'HAND1' %in% unlist(CTS)
#[1] TRUE

#### extract CTS scores for each biomodule candidate in the following steps: ####
## first to record the max MCI for the n.state.candidate 
maxMCI <- maxMCI[names(CTS)[1:n.state.candidate]]
maxMCI
#       9       10        7        1        4        8        5 
# 5.779724 5.276104 4.649539 4.151014 3.751412 3.640540 3.583947 

## then applendix the 2nd highest MCI score (if existing) for the states with exactly 2 bars
if(length(whoistop2nd)>0) maxMCI <- c(maxMCI, getNextMaxStats(membersL[['MCI']], idL=maxMCIms[['idx']], whoistop2nd))
names(maxMCI)
#[1] "9"  "10" "7"  "1"  "4"  "8"  "5"  "9"  "10" "7"  "1"  "8"  "5" 

# ## applendix the 2nd highest MCI score (if existing) for the states with exactly 3 bars
# if(length(whoistop3rd)>0) maxMCI <- c(maxMCI, getNextMaxStats(membersL[['MCI']], idL=maxMCIms[['idx']], whoistop3rd))
# names(maxMCI)
# # [1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2"   "Day1"   "Day2.5" "Day1.5" 
# 
# ## then applendix the 3rd highest MCI score (if existing) for the states with exactly 3 bars
# if(length(whoistop3rd)>0) maxMCI <- c(maxMCI, getNextMaxStats(membersL[['MCI']], idL=maxMCIms[['idx']], whoistop3rd, which.next=3))
# maxMCI
# #   Day2.5      Day2    Day1.5      Day1      Day2      Day1    Day2.5    Day1.5    Day2.5    Day1.5 
# # 4.2294562 3.9979599 3.6935297 3.6688229 1.9068669 0.7693816 3.5453475 3.0876365 3.4842829 1.0100347 


## to ensure the same order between maxMCI  and CTS
all(names(CTS) == names(maxMCI)) # TRUE


save(CTS,  maxMCI, whoistop2nd, file=paste0("CTS_maxMCI.RData")) # !!!!!!!!!!!!!!

## estimate significance NOT repeat (takes a while)
C = 1000

simuMCI <- list()
for(i in 1:length(CTS) ){
  n <- length(CTS[[i]])
  simuMCI[[i]] <- simulationMCI(n, samplesL, logmat,  B=C, fun="BioTIP")
}
names(simuMCI) <- names(CTS)


## the sig CTS candidates were c(1:4,7,9)
pdf(file=paste0("MCI_GenePermutation_1000_CTS_timepoint", cut.preselect, "_fdr",cut.fdr,"_minsize",cut.minsize,".pdf"),
    width=10)
P <- maxMCI * 0
par(mfrow=c(2,4))
for(i in 1:length(CTS) ){
  n <- length(CTS[[i]])
  P[i] <- length(which(maxMCI[i] <= simuMCI[[i]][names(CTS)[i],]))/C
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las=2, ylim=c(0,6),
                    main=paste(names(CTS)[i], n, "genes",
                               "\n","vs. ",C, "times of gene-permutation"), which2point=names(CTS)[i])
}
dev.off()
P
# 9 10  7  1  4  8  5  9 10  7  1  8  5 
# 0  0  0  0  0  0  0  0  0  0  0  0  0 
(sig <- which(P<0.01) ) 
#9 10  7  1  4  8  5  9 10  7  1  8  5 
#1  2  3  4  5  6  7  8  9 10 11 12 13 

save(simuMCI,P, file=paste0("GenePermutation_",C,"CTS_timepoint", cut.preselect, "_fdr",cut.fdr,"_minsize",cut.minsize,".RData"))


####### Verifying Tipping Point using IC* (or old Ic) score  #################
######  BioTIP score, permutating genes ##############
C= 1000
file=paste0("GenePermutation_",C,"CTS_timepoint", cut.preselect, "_fdr",cut.fdr,"_minsize",cut.minsize,".RData")
load(file)

load(file=paste0("CTS_maxMCI.RData"))

set.seed(2021)
IC_new <- simuBioTIP_g <- list()  
for(i in 1:length(sig)){
  n <- length(CTS[sig][[i]])
  IC_new[[i]] <- getIc(logmat, samplesL, CTS[sig][[i]], fun="BioTIP", shrink=TRUE)

  simuBioTIP_g[[i]]  <- simulation_Ic(n, samplesL, logmat, B=C,
                                 fun="BioTIP", shrink=TRUE)
}
names(simuBioTIP_g) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
names(IC_new) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')


names(CTS)[sig]
#[1] "9"  "10" "7"  "1"  "4"  "8"  "5"  "9"  "10" "7"  "1"  "8"  "5" 

pdf(file="BioTIP_vsSimulation_13CTS.pdf", width = 10)
par(mfrow=c(3,4))
for(i in 1:length(sig)){
  plot_Ic_Simulation(IC_new[[i]], simuBioTIP_g[[i]], las = 2, ylab="Ic*",
                     main= names(IC_new)[i],
                     fun="boxplot", which2point= names(CTS)[sig][i],
                     ylim=c(0,0.4)) # matplot
 }
dev.off()

P.g_simu <- array()
par(mfrow=c(3,4))
for(i in 1:length(sig)){
   P.g_simu[i] <- plot_SS_Simulation(IC_new[[i]], simuBioTIP_g[[i]], 
                                    main = paste("Delta Ic*"), 
                                    ylab=NULL)
}


P.g_simu
#[1] 0.000 0.000 0.011 0.165 0.423 0.001 0.000 0.000 0.179 0.831 0.558 0.252 0.257
(sig <- which(P.g_simu<0.05))
#[1] 1 2 3 6 7 8

save(IC_new, simuBioTIP_g, P.g_simu,  file="simuBioTIP_g.RData", compress=TRUE)
unlist(lapply(IC_new, max))
#   9_18g    10_19g     7_10g     1_19g     4_23g     8_16g     5_27g     9_23g 
# 0.3435304 0.3681858 0.3023271 0.1814797 0.2615350 0.2583500 0.2907017 0.2984304 
#   10_14g     7_23g     1_14g     8_38g     5_29g 
# 0.1743083 0.1917179 0.2076771 0.1939648 0.1699749 

  
# ######  old IC score, shulffing genes ##############
# set.seed(2021)
# IC_old <- simuIc_g <- list()  
# for(i in 1:length(sig)){
#   n <- length(CTS[sig][[i]])
#   IC_old[[i]] <- getIc(logmat, samplesL, CTS[sig][[i]], fun="cor", shrink=TRUE)
#   
#   simuIc_g[[i]]  <- simulation_Ic(n, samplesL, logmat, B=C,
#                                       fun="cor", shrink=TRUE)
# }
# names(simuIc_g) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
# names(IC_old) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
# 
# save(IC_old, simuIc_g,   file="simuIC_g.RData", compress=TRUE)
# 
# names(CTS)[sig]
# #[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2.5" "Day2.5"  
# 
# pdf(file="IC_vsSimulation_6CTS.pdf", height=10)
# par(mfrow=c(3,4))
# for(i in 1:length(IC_old)){
#   plot_Ic_Simulation(IC_old[[i]], simuIc_g[[i]], las = 2, ylab="old Ic",
#                      main= names(IC_new)[i],
#                      fun="boxplot", which2point= names(CTS)[sig][i]) # matplot
#   plot_SS_Simulation(IC_old[[i]], simuIc_g[[i]], 
#                      main = paste("Delta Ic"), 
#                      ylab=NULL)
#   
# }
# dev.off()

####### Verifying Tipping Point and using IC* (or old Ic) score  #################
######  BioTIP score, Shuffling label ##############
C= 1000
lengths(samplesL)
#  1  10   2   3   4   5   7   8   9 
# 77 107 161  86  86  80 173  95  49 


set.seed(2021)
simuBioTIP_s <- list()  

for(i in 1:length(sig)){
 # IC_new[[i]] <- getIc(logmat, samplesL, CTS[sig][[i]], fun="BioTIP", shrink=TRUE)
  simuBioTIP_s[[i]] <- matrix(nrow=length(samplesL), ncol=C)
 # rownames(simuBioTIP_s[[i]]) = names(CTS)[sig][i]
  rownames(simuBioTIP_s[[i]]) = names(samplesL)
  for(j in 1:length(samplesL)) { 
    #ns <- length(samplesL[names(CTS[sig])][[i]])  # for this state only
    ns <- length(samplesL[[j]])  # for each state rewpectively 
    simuBioTIP_s[[i]][j,]  <- simulation_Ic_sample(logmat, ns, Ic=IC_new[[i]],
                                             genes=CTS[sig][[i]], B=C,
                                             fun="BioTIP", shrink=TRUE)
  }
}
names(simuBioTIP_s) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
#names(IC_new) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')


names(CTS)[sig]
#[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2.5" "Day2.5"  

P.s_simu <- array(dim=length(sig))

pdf(file=paste0("BioTIP_vs_Shufflinglabel_",length(sig), "CTS.pdf"), width=10)
par(mfrow=c(3,4))
for(i in 1:length(sig)){
  plot_Ic_Simulation(IC_new[sig][[i]], simuBioTIP_s[[i]], las = 2, ylab="Ic*",
                     main= names(IC_new)[i],
                     fun="boxplot", which2point= names(CTS)[sig][i]) # matplot
  P.s_simu[i] <- plot_SS_Simulation(IC_new[sig][[i]], simuBioTIP_s[[i]], 
                     main = paste("Delta Ic*"), 
                     ylab=NULL)
  
}
dev.off()

P.s_simu
#[1] 0.001 0.000 0.093 0.005 0.011 0.000

save(simuBioTIP_s, P.s_simu,  file="simuBioTIP_s.RData", compress=TRUE)


CTS[sig]
# $`9`
# [1] "DKK1"   "DLL3"   "FGF8"   "FOXH1"  "FZD4"   "FZD7"   "GATA6"  "GSC"    "HRT2"  
# [10] "KDR"    "KIT"    "MESP1"  "MESP2"  "MIXL1"  "NKX2.5" "T"      "TBX5"   "SHH"   
# 
# $`10`
# [1] "BMP2"  "DLL1"  "DLL3"  "EVX1"  "FGF12" "FGF8"  "FOXC1" "FZD7"  "GATA4" "HAND2" "HRT2" 
# [12] "ISL1"  "MESP1" "MESP2" "MIXL1" "SFRP1" "TBX2"  "WNT4"  "WNT5A"
# 
# $`7`
# [1] "BMP4"  "DLL1"  "DLL3"  "FOXC1" "FZD2"  "HAND1" "HRT2"  "MESP2" "TBX1"  "WNT5A"
# 
# $`8`
# [1] "BAMBI"  "DKK1"   "EVX1"   "HEY1"   "HNFA4"  "HRT2"   "MESP1"  "MESP2"  "MSX1"  
# [10] "MSX2"   "MYL4"   "PDGFRA" "PDGFRB" "SERCA"  "T"      "WNT4"  
# 
# $`5`
# [1] "BMP2"  "CXCR4" "DKK1"  "DLL3"  "EVX1"  "FGF10" "FGF12" "FGF8"  "FZD2"  "GAS1"  "GATA4"
# [12] "GATA6" "GSC"   "HNFA4" "INHBA" "KIT"   "MIXL1" "MSX1"  "PTX2"  "SFRP1" "T"     "SOX17"
# [23] "WNT3A" "WNT5B" "TBX1"  "MYL3"  "TNNT2"
# 
# $`9`
# [1] "ACVR1B"  "ACVR2A"  "BMP2"    "BMP4"    "BMPR1A"  "EMILIN2" "FGF12"   "FSTL1"  
# [9] "FZD2"    "FZD6"    "GATA4"   "LTBP1"   "MSX1"    "MYOCD"   "NUMB"    "PARD3"  
# [17] "PDGFA"   "PDGFB"   "PDGFRB"  "RPL35A"  "VEGFA"   "SIRPA"   "WNT5A"  


######  old IC score, shulffing labels ##############
# set.seed(2021)
# IC_old <- simuIc_s <- list() 
# 
# for(i in 1:length(sig)){
#   IC_old[[i]] <- getIc(logmat, samplesL, CTS[sig][[i]], fun="cor", shrink=TRUE)
#   simuIc_s[[i]] <- matrix(nrow=length(samplesL), ncol=C)
#   rownames(simuIc_s[[i]]) = names(samplesL)
#   for(j in 1:length(samplesL)) { 
#     ns <- length(samplesL[[j]])  # for each state rewpectively 
#     simuIc_s[[i]][j,]  <- simulation_Ic_sample(logmat, ns, Ic=IC_new[[i]],
#                                                  genes=CTS[sig][[i]], B=C,
#                                                  fun="cor", shrink=TRUE)
#  }
# }
# names(simuIc_s) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
# names(IC_old) <- paste0(names(CTS)[sig],'_',lengths(CTS)[sig],'g')
# 
# save(IC_old, simuIc_s,   file="simuIc_s.RData", compress=TRUE)
# 
# names(CTS)[sig]
# #[1] "Day2.5" "Day2"   "Day1.5" "Day1"   "Day2.5" "Day2.5"  
# pdf(file="IC_vs_Shufflinglabel_6CTS.pdf", height=10)
# par(mfrow=c(3,4))
# for(i in 1:length(IC_old)){
#   plot_Ic_Simulation(IC_old[[i]], simuIc_s[[i]], las = 2, ylab="old Ic",
#                      main= names(IC_new)[i],
#                      fun="boxplot", which2point= names(CTS)[sig][i]) # matplot
#   plot_SS_Simulation(IC_old[[i]], simuIc_s[[i]], 
#                      main = paste("Delta Ic"), 
#                      ylab=NULL)
#   
# }
# dev.off()
# 
# 

