setwd("E:/Git_Holly/scRNAseq_examples/result/")
source(file='E:/Git_Holly/scRNAseq_examples/code/Stability_score.R')

## namually record teh Ic predictions

subfolds <- c('hESC_Bargaje2017', 'lung_Treutlein2014', 
              'gastrulationE8.25_Pijuan-Sala2019',
              'gastrulationE8.25_Ibarra-Soria2018',
              'EB_Zhao2019','simulated_EMT')

DNB <- list(hESC='C9', # if using the highest score in the system without evaluating significance 
                lung= 'C2',
                E8.25_2019 = 'C15',  
                E8.25_2018 = 'endothelial.b',
                EB = 'C5',
                EMT = 'TC')  

Ic <- list(hESC='C9',
           lung='C5',
           E8.25_2019 = c('C13','C15'),
           E8.25_2018 = c('endothelial.b', 'endothelial.d'),
           EB = c('C10','C17','C4'),
           EMT = 'TC')
Dbs <- names(Ic)

## according to the results of stability testing in each dataset
CT.ref <- list(hESC=c('C9','C10'),
               lung='C2',
               E8.25_2019 = c('C13','C15','C6','C7'),
               E8.25_2018 = c('endothelial.b', 'cardiac.a'),
               EB = c('C11','C6','C4'),
               EMT = c('TC','I1'))
names(CT.ref) <- Dbs

## reead the BioTIP's identification ---------------
cluster.methods <- c('C_Consensus_929cells', 'C_Leiden_0.4', 
                     'C_SNNGraph', 'subcelltype',
                     'C_SNNGraph_k5', 'C_cell_type')
BioTIP <- list()
for(i in 1:length(subfolds)){
  if(i==3) {
   # load(file= paste0(subfolds[i], '/', cluster.methods[i],'/', 'CTS.RData'))
   # BioTIP[[i]] <- paste0('C',unique(names(CTS)))
    BioTIP[[i]] = c("C15", "C6", "C13") # C16 was excluded in the stability testing uisng only subset of datasets
    ## therefore S16 was a negative control in the stability test
  } else {
    load(file= paste0(subfolds[i], '/', cluster.methods[i],'/', 'BioTIP.res.RData'))
    CT <- names(res$CTS.candidate)[which(res$significant)]
    if(!grepl('C', CT[1]) & !grepl('.', CT[1], fixed=T)) CT <- paste0('C', CT)
    BioTIP[[i]] <- unique(CT) 
  }  
}
names(BioTIP) <- Dbs

## reead the QuanTC's identification ---------------
QuanTC <- list()
for(i in 1:length(subfolds)){ 
  if(i==1)  QuanTC[[i]] <- c('C9', 'C10')
  if(i==2)  QuanTC[[i]] <- c('C2','C4') # both E16.5 and E14.5 cells were identified but E16.5 is the highest, 
  if(i==3)  QuanTC[[i]] <- c('C15','C6', 'C13','C3') # using the better outputs of QuanTC ran on (k=6)
  if(i==4)  QuanTC[[i]]=NA # no QuanTC outputs for E8.25 2018 data
  if(i==5)  QuanTC[[i]] <- c('C1','C4', 'C2', 'C6', 'C5', 'C8', 'C3', 'C7')
#     {
#     load(file= paste0(subfolds[i], '/QuanTC_run/CTS.RData'))
#     CT <- names(CTS)
# #    if(!grepl('C', CT) & !grepl('.', CT, fixed=T)) CT <- paste0('C', CT)
#     if(!grepl('-', CT, fixed=TRUE)) CT <- lapply(CT, function(x) unlist(strsplit(x, '-'))) %>% unlist() %>% unique()
#     QuanTC[[i]] <- unique(CT) 
#   }
  if(i==6)  QuanTC[[i]] <- c('TC','I1')
  
}  
names(QuanTC) <- Dbs

## calculate the jaccard scores --------------------------
jaccard.CT <- matrix(nrow=4, ncol=length(Dbs))
row.names(jaccard.CT) <- c('DNB','IC', 'QuanTC','BioTIP')
colnames(jaccard.CT) <- Dbs
for(i in 1:length(Dbs)){
  jaccard.CT[1,i] <- jaccard.sim(CT.ref[[i]], DNB[[i]]) 
  jaccard.CT[2,i] <- jaccard.sim(CT.ref[[i]], Ic[[i]]) 
  jaccard.CT[3,i] <- jaccard.sim(CT.ref[[i]], QuanTC[[i]]) 
  jaccard.CT[4,i] <- jaccard.sim(CT.ref[[i]], BioTIP[[i]]) 
}

### bar plot Jarcarrd scores ---------------
df <- data.frame( method = rep(rownames(jaccard.CT), ncol(jaccard.CT)),
                  Jaccard = as.numeric(jaccard.CT),
                  dataset = lapply(colnames(jaccard.CT), 
                                   function(x) rep(x, nrow(jaccard.CT))) %>% unlist())
df$method = factor(df$method, levels=c('BioTIP','QuanTC','IC','DNB'))
df$dataset = factor(df$dataset, levels= rev(Dbs))

p1 <- ggplot(data=df, aes(x=dataset, y=Jaccard, fill=method, width=.6)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #scale_fill_brewer(palette="Paired")+
  scale_fill_manual( values = c("DNB"="green4", "IC" = "sky blue", "QuanTC" = "blue", "BioTIP" = "orange")) +
  theme_minimal()  + coord_flip()  

#########################################################
## extract the BioTIP's F1 scores, comapring to GS or PGS per dataset --------------------------
cluster.to.calculate.F1 <- c('Norm.F1.PS', 'Norm.F1.CTS', 'Norm.F1.eHEP', 
                             'endothelial.b', 'Norm.F1.KLF1',
                             'Norm.F1.TC')
names(cluster.to.calculate.F1) <- Dbs

F1.norm <- matrix(nrow=2, ncol=length(Dbs))
row.names(F1.norm) <- c('QuanTC','BioTIP')
colnames(F1.norm) <- Dbs

for(i in 1:length(subfolds)){
  if(i==4) {tb <- read.table(file= paste0(subfolds[i], '/stability/jaccard.CT_F1.CTS.txt')) 
    F1 = c(tb$F1.End.ctl,   tb$F1.End,
         tb$F1.Cardiac.c.ctl,    tb$F1.Cardiac.c)
    Norm.F1 = Normalize.F1(F1)
    x <- length(F1)/4
    tb$Norm.F1.End <- Norm.F1[(x*1+1):(x*2)]
    F1.norm['BioTIP',i] <- mean(tb[, 'Norm.F1.End'])
    F1.norm['QuanTC',i] <- 0  
  } else {
      tb <- read.table(file= paste0(subfolds[i], '/jaccard.CT_F1.CTS.txt'))  
      if(i==3) rownames(tb)[which(rownames(tb)=='QuanTC.k6_run')] = 'QuanTC_run'
      F1.norm['BioTIP',i] <- tb[cluster.methods[i], cluster.to.calculate.F1[i]]
      F1.norm['QuanTC',i] <- tb['QuanTC_run', cluster.to.calculate.F1[i]]  
   }
}

F1.norm
#       hESC  lung E8.25_2019 E8.25_2018    EB       EMT
# QuanTC 0.651 0.322      0.307   0.000000 0.247 0.6154848
# BioTIP 1.000 0.952      0.485   0.753395 0.937 0.9833469

## extract the DNB's F1 scores --------------------------
# for which db3 (E8.25 2019), db5 (EB) data can not detect the cluster of interest with the highest DNB score
F1.norm <- rbind(F1.norm, 'DNB'=rep(0, ncol(F1.norm)))
#F1.norm['DNB', 'E8.25_2018'] <- 1 ## we used the C5_33g as gold standard
F1.norm['DNB', 'E8.25_2018'] <- F1.norm['BioTIP', 'E8.25_2018'] ## we used the C5_33g as gold standard


DNB.subfold.ID <-c(1,2,6)
cluster.to.calculate.F1[DNB.subfold.ID] <- c('Norm.F1.PS.DNB', 'Norm.F1.DNB','Norm.F1.TC')
for(i in DNB.subfold.ID){
  tb <- read.table(file= paste0(subfolds[i], '/jaccard.CT_F1.CTS.txt'))  
  F1.norm['DNB',i] <- tb[cluster.methods[i], cluster.to.calculate.F1[i]]
}
  

### bar plot normalized F1 scores ---------------
df2 <- data.frame( method = rep(rownames(F1.norm), ncol(F1.norm)),
                  F1.norm = as.numeric(F1.norm),
                  dataset = lapply(colnames(F1.norm), 
                                   function(x) rep(x, 3)) %>% unlist())
df2$method = factor(df2$method, levels=c('BioTIP','QuanTC','DNB'))
df2$dataset = factor(df2$dataset, levels= rev(Dbs))

p2 <- ggplot(data=df2, aes(x=dataset, y=F1.norm, fill=method, width=.6)) +
    geom_bar(stat="identity", position=position_dodge())+
   # scale_fill_brewer(palette="Paired")+
  scale_fill_manual( values = c("DNB"="green4", "QuanTC" = "blue", "BioTIP" = "orange")) +
    theme_minimal()  + coord_flip()  

pdf(file=paste0('Jarcard_CTS.F1_6dbs.pdf'), height=2)
gridExtra::grid.arrange(
  p1, p2, left = "CT detection", right = "CTS identification"
  ,ncol=2)
dev.off()
