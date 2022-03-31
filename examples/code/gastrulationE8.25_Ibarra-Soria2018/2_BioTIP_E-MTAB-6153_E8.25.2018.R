#### README ############
## There are 3 sections in this code
## section 0) plot TSNE for 1531 cells that we analyze
## Section 1) apply BioTIP to predefined subcelltypes of 11k cells             
## section 2) stability of BioTIP's identification
## Section 3) estimating the robustness using the downsized experiment

setwd('F:/projects/scRNA/results/IbarraSoria2018_MouseE8.25/E8.25.2018_robustness/')

# the following funciton has been updated on github, but before updating the package to Bioconductor, you need to call it locally to test
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/ForJennifer/optimize.sd_selection.R') 
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/BioTIP.wrap.R')
source(file='F:/projects/scRNA/source/GSE87038_gastrulation/GSE87038_git_copy/Stability_score.R')

  
# normalized data were downloaded and reanalyzed
# refer to GSE87038_E8.25_normalize.R
library(BioTIP)
library(dplyr)
library(SingleCellExperiment)

library(ggplot2)
library(gridExtra)  # required by grid.arrange()
library(ggpubr)  # required by p-values on bar

library(MouseGastrulationData)
head(EmbryoCelltypeColours)
length(EmbryoCelltypeColours) # 37


###################################################
## load the R object ##
## focusing on 131 cells along the AT2 trajectory
###################################################
load('sce_16subtype.RData')
sce
# class: SingleCellExperiment 
# dim: 4000 11039 
# metadata(0):
#   assays(3): counts normcounts logcounts
# rownames(20483): Xkr4 Gm1992 ... Csf2ra Gm21060
# rowData names(5): GENEID SYMBOL SEQNAME HVGs.10 HVGs.20
# colnames(11039): extraembryonicMesoderm_1 extraembryonicMesoderm_2 ...
# cardiac.c_281 cardiac.c_282
# colData names(3): label celltype subcelltype
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):
  
# color for plotting
MyEmbryoCelltypeColours <- EmbryoCelltypeColours[1:length(unique(sce$subcelltype))]
names(MyEmbryoCelltypeColours) <- levels(sce$subcelltype)


# logmat <- as.matrix(assays(sce)$logcounts)
# M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
# save(M, file="ShrinkM.RData", compress=TRUE)
# dim(M) #4000 4000
# rm(logmat)
load(file="ShrinkM.RData")
dim(M) #4000 4000


####################################################################
## section 1) running BioTIP on different predefiend cellsubtypes ##
####################################################################
{
  ######### 1.1) setting parameters,   ######################
  localHVG.preselect.cut = 0.1 # A positive numeric value setting how many propotion of global HVG to be selected per cell cluster
  getNetwork.cut.fdr = 0.05  # to construct RW network and extract co-expressed gene moduels
  getTopMCI.gene.minsize = 30  # min number of genes in an identified CTS
  getTopMCI.n.states = 9  # A number setting the number of states to check the significance of DNB score (i.e. MCI) 
  # This parameter can be enlarge when inputting more clusters, usually the expected number of transition states +1
  MCIbottom = 2  # A number setting the threshold to prioritize the initially selected CTS candidates.
  # In our experiment, a number between 2 and 4 generated expected resutls.
  
  ############### 1.2) applying BioTIP to all new clustering outputs  ##################
  colnames( colData(sce))
 # [1] "label"       "celltype"    "subcelltype"
  x <- grep("type", colnames(colData(sce)))
  
  
  set.seed(102020)
  for(i in 1:length(x)){
    samplesL <- split(rownames(colData(sce)), f = colData(sce)[x[i]])
    (tmp = lengths(samplesL))
    
    res <- BioTIP.wrap(sce, samplesL, subDir=colnames(colData(sce))[x[i]],
                       getTopMCI.n.states=getTopMCI.n.states, 
                       getTopMCI.gene.minsize=getTopMCI.gene.minsize, 
                       MCIbottom=MCIbottom,
                       local.HVG.optimize =TRUE, localHVG.preselect.cut=localHVG.preselect.cut,  
                       getNetwork.cut.fdr=getNetwork.cut.fdr, 
                       M=M,
                       permutation.method='gene',
                       verbose=TRUE, plot=TRUE)           
    save(res, file=paste0(colnames(colData(sce))[x[i]],'/BioTIP.res.RData'))  
    
  }
  
  ### control the max module size afterwards
  load(file='subcelltype/BioTIP.res.RData')
       x <- lengths(res$CTS.candidate)
       dropoff <- which(x> 100)
       if(length(dropoff)>0){
         res$CTS.candidate = res$CTS.candidate[-dropoff]
         res$topMCI = res$topMCI[-dropoff]
         res$BioTIP_scores = res$BioTIP_scores[-dropoff]
         res$significant = res$significant[-dropoff]
       }
   save(res, file='subcelltype/BioTIP.res.RData')
}


###########################################################
## Section 2) stability after downsized dataset          ##
## testing different getTopMCI.gene.minsize              ##
###########################################################
{
  load('sce_16subtype.RData')
  sce
  # class: SingleCellExperiment 
  # dim: 4000 11039 
  load(file="ShrinkM.RData")
  dim(M) #4000 4000

  colnames(colData(sce))
  # [1] "label"       "celltype"    "subcelltype"
  
  rowData(sce) <- NULL # to allow run sample()
  
  ######### 2.1) setting parameters, the same as section 1  ######################
  ## but try multiple getTopMCI.gene.minsize 
  getTopMCI.gene.minsize <- c(30,20,40)
  
  ######### 2.2) run BioTIP on downsized data (runs 1 day) ##################
  {
    n.scability <- 10
    downsize <- 0.95  
    n1 <- nrow(sce)
    n2 <- ncol(sce)
    
    #for(i in 1:length(MCIbottom)){
    i=1  
    set.seed(102020)
    ## get n.scability sets of predictions
    res <- list()  
    res[[1]] <- res[[2]] <- res[[3]] <- list()
    names(res) <- getTopMCI.gene.minsize
    
    for(flag in 1:n.scability){
      cat(flag, '\t')
      
      sce.dnsize <- sce[sample(n1, floor(n1*downsize)), sample(n2, floor(n2*downsize))]
      samplesL <- split(rownames(colData(sce.dnsize)), f = sce.dnsize$subcelltype)
      
      for(m in 1:length(getTopMCI.gene.minsize)){
        if(m>1){
          outputpath = paste0('stability/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
          load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
          logmat.local.HVG.testres = testres
          for(j in 1:length(testres)) samplesL[[j]] <- colnames(testres[[j]])
        } else logmat.local.HVG.testres = NULL
        this.run <- BioTIP.wrap(sce.dnsize, samplesL, subDir='stability', 
                                getTopMCI.n.states=getTopMCI.n.states, 
                                getTopMCI.gene.minsize=getTopMCI.gene.minsize[m],
                                #getTopMCI.gene.maxsize = 100,  # control the max module size because there are poulation over 1k cells
                                MCIbottom=MCIbottom[i],
                                localHVG.preselect.cut=localHVG.preselect.cut, 
                                logmat.local.HVG.testres = logmat.local.HVG.testres,
                                getNetwork.cut.fdr=getNetwork.cut.fdr, 
                                M=M,
                                permutation.method='gene',
                                verbose=TRUE, plot=FALSE)   
        res[[m]][[flag]] <- this.run
      }  # end loop m
    }   # end loop n.scability
    
    res.bk <- res
    for(m in 1:length(getTopMCI.gene.minsize)){
      res <- res.bk[[m]]
      save(res, flag, counts, file=paste0('stability/BioTIP_res_10runs_MCIbottom',MCIbottom[i],
                                          '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    }
  }

  
 ## Alternatively, control the max module size after first running without this control
  for(m in 1:length(getTopMCI.gene.minsize)){
    load(file=paste0('stability/BioTIP_res_10runs_MCIbottom2_Modulesize',getTopMCI.gene.minsize[m],'.RData'))

    for(i in 1:length(res)){
      x <- lengths(res[[i]]$CTS.candidate)
      dropoff <- which(x> 100)
      if(length(dropoff)>0){
        res[[i]]$CTS.candidate = res[[i]]$CTS.candidate[-dropoff]
        res[[i]]$topMCI = res[[i]]$topMCI[-dropoff]
        res[[i]]$BioTIP_scores = res[[i]]$BioTIP_scores[-dropoff]
        res[[i]]$significant = res[[i]]$significant[-dropoff]
      }
    }
    save(res, file=paste0('stability/BioTIP_res_10runs_MCIbottom2_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    
  }
   
  ## 2.3) evaluate the stability of CT detection
  ## using the robustly detected C2 as reference
  ###########################################################
  df.CT =NULL
  for(m in 1:length(getTopMCI.gene.minsize)){
    CT.list <- list()
    Freq <- matrix(0, nrow=nlevels(sce$subcelltype), ncol=length(MCIbottom))#!!!
    rownames(Freq) <- levels(sce$subcelltype) #!!!
    colnames(Freq) <- MCIbottom
    
    #  for(i in 1:length(MCIbottom)){
    i=1
    load(file=paste0('stability/BioTIP_res_10runs_MCIbottom',MCIbottom[i],
                     '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    ## count which cluster has been identified repeatily ###########
    CT <- list()
    # for(j in 1:n.stability)  
    for(j in 1:length(res)){ # in case n.stability>length(res), i.e. some run has no isgnificant prediction
      CT[[j]] <- res[[j]]$significant[which(res[[j]]$significant)] %>% names() %>% unique()
    }
    tmp <- table(unlist(CT))
    Freq[names(tmp),i] <- tmp
    #Freq <- cbind(Freq, tmp)
    CT.list[[m]] <- CT
    #  } 
    
    # Freq <- Freq/length(res)
    names(CT.list) <- MCIbottom 
    
    df <- data.frame(Frequency = as.numeric(Freq),
                     clusterID = rep(rownames(Freq), ncol(Freq)),
                     MCIbottom = lapply(as.character(MCIbottom ), function(x) rep(x,nrow(Freq))) %>% unlist()
    )
    head(df)
    #   Frequency clusterID MCIbottom
    # 1    10       C13         2
    # 2     1        C3         2
    # 3     5        C5         2
    # 4     2        C7         2
    # 5     3        C8         2
    # 6    12        C6         2
    df$clusterID <- factor(df$clusterID, levels =levels(sce$subcelltype)) #!!!
    df$moduleSize <- rep(getTopMCI.gene.minsize[m], nrow(df))
    
    df.CT <- rbind(df.CT, df)  
    save(CT.list, file= paste0('stability/CT.list_Modulesize',getTopMCI.gene.minsize[m],'.Rdata'))
  }
  
  df.CT$Frequency = df.CT$Frequency/n.scability
  save(df.CT, 
       file= 'stability/CT.detection_variable.Modulesize.Rdata')
  
  
  ## 2.4) plot for variable Modulesize, when MCIbottom==2         ##
  #############################################################
  load(file= 'stability/CT.detection_variable.Modulesize.Rdata')
  df.CT$moduleSize <- as.character(df.CT$moduleSize )
  df.CT$clusterID <- factor(df.CT$clusterID,
                            levels = c('mesodermProgenitors', 'presomiticMesoderm.b','somiticMesoderm',
                                       'presomiticMesoderm.a', 'mixedMesoderm.a','mixedMesoderm.b',
                                       'pharyngealMesoderm', 'extraembryonicMesoderm',
                                       'cardiac.a', 'cardiac.b', 'cardiac.c',
                                       'endothelial.a','endothelial.b','endothelial.c','endothelial.d',
                                       'blood'))
  ggplot(data=subset(df.CT, MCIbottom==2),
         aes(x=clusterID, y=Frequency, group=moduleSize)) +
    geom_line(aes(color=moduleSize), size=1.5)+
    geom_point(aes(color=moduleSize)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_color_manual(values=c( "orange","green","blue")) 
  dev.copy2pdf(file='stability/CT.detection_variable.Modulesize_MCIbottom2.pdf')
  
  
  ## 2.5) evaluate the stability of CTS identification
  ## for the three detected subtypes, respectively
  ## between the downsized prediction and fullsized prediction
  ###########################################################
  load('stability/CT.list_Modulesize30.Rdata')
  names(CT.list)
  #[1] "2"  
  
  load('subcelltype/BioTIP.res.RData')
  res$CTS.candidate[which(res$significant)] %>% names()
  #[1] "endothelial.b" "cardiac.a"         
  (End.ref = res$CTS.candidate[1])
  # $endothelial.b
  # [1] "Etv2"     "Tal1"     "Hhex"     "Ctla2a"   "Gngt2"    "Sox17"   
  # [7] "Ifitm3"   "Cd34"     "Icam2"    "Cdh5"     "Ebf1"     "Rhoj"    
  # [13] "Igf1"     "Grrp1"    "Abi3"     "Gata2"    "Gimap1"   "Fxyd5"   
  # [19] "Nrp2"     "Gcc2"     "Gab1"     "Fev"      "Fam89a"   "Mageh1"  
  # [25] "Slc66a2"  "Lhfp"     "Dusp4"    "Npr1"     "C1ql2"    "Tbc1d22b"
  # [31] "Col13a1"  "Cd40"     "Tek" 
  (Cardiac.a.ref = res$CTS.candidate[2])
   
  
  MCIbottom  # 2
  i=1
  
  res.End <- get.downsized.F1.score(End.ref, getTopMCI.gene.minsize, MCIbottom, i=1)
  save(res.End, file='stability/F1.CTS.endothelial.b_MCIbottom2.RData')
  
  res.Cardiac.a <- get.downsized.F1.score(Cardiac.a.ref, getTopMCI.gene.minsize, MCIbottom, i=1)
  save(res.Cardiac.a, file='stability/F1.CTS.eCardiac.a_MCIbottom2.RData')
  
  
  lapply(res.End$F1, mean) %>% unlist()
  #       30        20        40 
  # 0.3628396 0.3789869 0.2793133 
  
  df <- data.frame(F1 = c(unlist(res.End$F1.ctl),unlist(res.End$F1),
                          unlist(res.Cardiac.a$F1.ctl),unlist(res.Cardiac.a$F1) ),
                   minModuleSize = c(rep(getTopMCI.gene.minsize, lengths(res.End$F1.ctl)),
                                     rep(getTopMCI.gene.minsize, lengths(res.End$F1)),
                                     rep(getTopMCI.gene.minsize, lengths(res.Cardiac.a$F1.ctl)),
                                     rep(getTopMCI.gene.minsize, lengths(res.Cardiac.a$F1)) ),
                   CT = c(rep('End', sum(lengths(res.End$F1.ctl))+sum(lengths(res.End$F1))),
                          rep('Cardiac.a', sum(lengths(res.Cardiac.a$F1.ctl))+sum(lengths(res.Cardiac.a$F1))) ),
                   CTS = c(rep('ctl', sum(lengths(res.End$F1.ctl))), rep('CTS', sum(lengths(res.End$F1))),
                           rep('ctl', sum(lengths(res.Cardiac.a$F1.ctl))), rep('CTS', sum(lengths(res.Cardiac.a$F1))) ) 
  )
  df$minModuleSize <- as.character(df$minModuleSize) 
  df$CTS <- factor(df$CTS, levels = c('ctl','CTS')) 
  df$color.lab <- paste0(df$minModuleSize,df$CTS)
  df$color.lab <- factor(df$color.lab, levels=c('20ctl','20CTS','30ctl','30CTS','40ctl','40CTS')) #, '40ctl','40CTS'))
  df$CT <- factor(df$CT, levels =c('End','Cardiac.a'))
  
  df$Norm.F1 = Normalize.F1(df$F1)
  
  # To ensure that easy and difficult CT state have equal influence on the final score, 
  # we normalized the scores at each CT state across the different MinModuleSize setting.
  
  library(ggbeeswarm)
  library(grid)
  library(gridExtra)
  
  # plot w C6, used 
  {
    p1 <- ggplot(data=df, aes(x=CT, y=F1, color=color.lab)) +
      geom_violin(position = position_dodge(width = 0.5)) +
      geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
      xlab('predicted for cluster') + ylab('F1 score for the CTS') +
      scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue", "grey","red")) 
    
    p2 <- ggplot(data=df, aes(x=CT, y=Norm.F1, color=color.lab)) +
      geom_violin(position = position_dodge(width = 0.5)) +
      geom_quasirandom(dodge.width = 0.5, varwidth = TRUE) + 
      xlab('predicted for cluster') + ylab('Normalized F1 score for the CTS') +
      scale_color_manual(values=c("grey", "orange","grey", "green","grey","blue", "grey","red")) 
    
    pdf(file=paste0('stability/CTS.F1_variable.ModuleSize_GS_2CTs.pdf'))
    gridExtra::grid.arrange(
      p1, p2, top = "CTS stability", bottom = "different ModuleSize settings"
      ,ncol=1)
    dev.off()
  }
  
  dim(sce)  # 4000 11039
  ## prsent the statistics using proxy gold standard here an in the figure
  # ns: p > 0.05
  # *: p <= 0.05
  # **: p <= 0.01
  # ***: p <= 0.001
  # ****: p <= 0.0001
  
  tmp = subset(df, minModuleSize==20 & CT=='End')
  t.test(tmp$F1[which(tmp$CTS=='CTS')], tmp$F1[which(tmp$CTS=='ctl')])
  #t = 7.0921, df = 9.0107, p-value = 5.682e-05   ****
  tmp = subset(df, minModuleSize==30 & CT=='End')
  t.test(tmp$F1[which(tmp$CTS=='CTS')], tmp$F1[which(tmp$CTS=='ctl')])
  #t = 5.435, df = 9.0238, p-value = 0.0004099  ***
  tmp = subset(df, minModuleSize==40 & CT=='End')
  t.test(tmp$F1[which(tmp$CTS=='CTS')], tmp$F1[which(tmp$CTS=='ctl')])
  #t = 3.3908, df = 9.0439, p-value = 0.007933  **
  
  tmp = subset(df, minModuleSize==20 & CT=='Cardiac.a')
  t.test(tmp$F1[which(tmp$CTS=='CTS')], tmp$F1[which(tmp$CTS=='ctl')])
  #t = 4.2209, df = 9.0244, p-value = 0.002223  **
  tmp = subset(df, minModuleSize==30 & CT=='Cardiac.a')
  t.test(tmp$F1[which(tmp$CTS=='CTS')], tmp$F1[which(tmp$CTS=='ctl')])
  #t = 2.917, df = 9.029, p-value = 0.01706   *
  tmp = subset(df, minModuleSize==40 & CT=='Cardiac.a')
  t.test(tmp$F1[which(tmp$CTS=='CTS')], tmp$F1[which(tmp$CTS=='ctl')])
  #t = 1.0232, df = 9.0067, p-value = 0.3329  ns
  
 
  
}


###################################################################################
## Section 3) estimating the robustness using the downsized experiment             ##
###################################################################################
{
  ## 3.1) calculate Jaccard-similarity score for identifing the 3 detected bufurcations 
  ## between each prediction of downsized data (n=10) and the full data (n=1; getTopMCI.gene.minsize=30; MCIbottom=2) --------------
  reference.full <- c("endothelial.b" ,"cardiac.a") 
  
  load(file=paste0('stability/CT.list_Modulesize30.RData'))
  jaccard.CT <- array()
  for(i in 1:n.scability){
    sub.CT <- CT.list[['2']][[i]]
    jaccard.CT[i] <- jaccard.sim(reference.full, sub.CT) 
  }
  round(jaccard.CT, 2)
  #[1] 0.17 0.33 0.50 0.25 0.25 0.67 0.67 0.67 0.33 0.20 
  
  ## 3.2) calculate F1 score 
  # 1 result from the full dataset ------------
  load('subcelltype/BioTIP.res.RData')
  CTSs <- res$CTS.candidate[which(res$significant)]
 
  # 10 results from the downsized datasets ------------
  load(file="stability/BioTIP_res_10runs_MCIbottom2_Modulesize30.RData")
  lengths(res)  #10
  F1 <- F1.ctl <- list()
  for(i in 1:length(CTSs)){
    F1.CTS <- F1.bk.CTS <- array()
    set1 <- CTSs[i]  
    for(j in 1:length(res)){
        set2 <- res[[j]]$CTS.candidate[which(res[[j]]$significant)]
        F1.CTS[j] <- F_score.CTS(set1, set2, weight=TRUE)$F1
        # negative control
        n <- lengths(set2)
        random.set2 <- lapply(n, function(x) sample(rownames(sce), x))
        F1.bk.CTS[j] <- F_score.CTS(set1, random.set2, weight=TRUE)$F1
    } 
    F1[[i]] <- F1.CTS
    F1.ctl[[i]] <- F1.bk.CTS
    }
  names(F1) <- names(F1.ctl) <- names(CTSs)  
  }

## which run identified Etv2?
HEP.marker <-  c('Etv2','Tal1','Rhoj','Rasip1','Lyl1','Plvap',)   

HEP.CTS <- NULL 
for(j in 1:length(res)){
  set2 <- res[[j]]$CTS.candidate[which(res[[j]]$significant)]
  if('endothelial.b' %in% names(set2)) tmp <- HEP.marker %in% unlist(set2[['endothelial.b']]) else tmp = FALSE
  HEP.CTS = rbind(HEP.CTS, tmp)
  }
colnames(HEP.CTS) <- HEP.marker
#       Etv2  Tal1  Rhoj
# tmp  TRUE  TRUE  TRUE
# tmp  TRUE  TRUE  TRUE
# tmp  TRUE  TRUE  TRUE
# tmp  TRUE  TRUE  TRUE
# tmp FALSE FALSE FALSE
# tmp FALSE FALSE FALSE
# tmp  TRUE  TRUE  TRUE
# tmp FALSE FALSE FALSE
# tmp  TRUE  TRUE  TRUE
# tmp FALSE FALSE FALSE

## reproting table #######################
tb <- data.frame(Jaccard.CT = round(jaccard.CT,3), 
                 F1.End.ctl = round(F1.ctl[[1]],3), F1.End = round(F1[[1]],3), 
                 F1.Cardiac.c.ctl = round(F1.ctl[[2]],3), F1.Cardiac.c = round(F1[[2]],3))
tb <- cbind(tb, HEP.CTS)
nrow(tb) #[1] 10

tb$detected.CT <- lapply(CT, toString) %>% unlist()

write.table(tb, file='stability/jaccard.CT_F1.CTS.txt',  sep='\t')


###############
## PLOT #######
###############

## generate bar plot, the first green bar is the average Jaccard.CT #######################
tb <- read.table('stability/jaccard.CT_F1.CTS.txt', header=T, sep='\t')
runs <- c(1:nrow(tb),11)
Jaccard.CT <- c(tb$Jaccard.CT, mean(tb$Jaccard.CT))
tmp <- data.frame(runs=runs,Jaccard.CT=Jaccard.CT )
  
p1 <- ggbarplot(tmp, "runs", "Jaccard.CT", orientation = "horiz") +
  geom_col(aes(fill = Jaccard.CT)) + 
  scale_fill_gradient2(low = "white",high = "green") +
  theme_bw(base_size=10)  +        # use a black-and-white theme with set font size
  theme(legend.position = "none") + coord_flip()   

## generate the bar plot for normalized F1 scores, when modulesize==30
df$aggregate.group <- paste(df$CT, df$CTS, df$minModuleSize,sep="_")
x <- split(df$Norm.F1, df$aggregate.group)
y <- lapply(x, mean) %>% unlist()
tmp <- data.frame(ave.Norm.F1=y,
                  ID = names(y),
                  CT = lapply(names(y), function(x) unlist(strsplit(x, "_"))[1]) %>% unlist(),
                  CTS = lapply(names(y), function(x) unlist(strsplit(x, "_"))[2]) %>% unlist(),
                  minModuleSize =lapply(names(y), function(x) unlist(strsplit(x, "_"))[3]) %>% unlist()
)

p2 <- ggbarplot(subset(tmp,  minModuleSize ==30), 
                "ID", "ave.Norm.F1", orientation = "horiz") +
  geom_col(aes(fill = ave.Norm.F1)) + 
  scale_fill_gradient2(low = "white",high = "blue") +
  theme_bw(base_size=10)  +        # use a black-and-white theme with set font size
  theme(legend.position = "none") + coord_flip()   


pdf(file='bar_robustness_fr.downSizedRuns.pdf', width=6, height=4)
gridExtra::grid.arrange(
  p1, p2#, p3  
  ,ncol=2)
dev.off()

