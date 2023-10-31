# This code uses the same methodology described in a previous benchmark (Saelens, Cannoodt et al. 2019) 
# to evaluate the stability of CTS identification using BioTIP.
# It calcualtes F score for two CTSs and
# summarizes the normalized stability score
# last update 1/17/2022
# Authos: Holly Yang (xyang2.at.uchicago.edu)
# last update: 01/17/2022

require(dplyr)

jaccard.sim = function(CTSx, CTSy) {
  CTSx <- unique(CTSx)
  CTSy <- unique(CTSy)
  y = length(intersect(CTSx, CTSy)) / length(union(CTSx, CTSy))
 # cat(length(intersect(CTSx, CTSy)) , '\n')
  return(y)
}

F_score.CTS <- function(Set1, Set2, weight=TRUE, exclude.unique.CT=TRUE, merge.cluser=FALSE){
  # Set1 is a named list of identified CTSs 
  # Set2 is another named list of identified CTSs
  # weight A logical value, default is TRUE, to leverage whether the Recover/Relevance is calcualted beween the same (or best-matched) clusters of two clustering outputs
  # exclude.unique.CT A logical value, default is TRUE, to compare the F1 score for commonly identifed CTs from two BioTIP identifications
  ## Otherwise, compare all identifications
  # merge.cluser A logical value, default is FALSE distincting multiple CTS identified for the same cluster ID, becasue BioTIP orders muterple CTSs for the same cluster. 
  ## However, this parameter will be set to FALSE when compare resutls of different methods (e.g., BioTIP vs QuanTC) because the order is nolonger comparable.
  if(any(is.na(names(Set1)))) stop('pls name you Set1 by its indeicated CT cluster ID')
  if(any(is.na(names(Set2)))) stop('pls name you Set2 by its indeicated CT cluster ID')
  
  if(exclude.unique.CT){
    common.CT <- intersect(names(Set1), names(Set2))
    if(length(common.CT)>0) {
      Set1 <- Set1[which(names(Set1) %in% common.CT)]
      Set2 <- Set2[which(names(Set2) %in% common.CT)]
    }
  }

  n1 <- length(Set1)
  n2 <- length(Set2)
  cluster1 <- unique(names(Set1))
  cluster2 <- unique(names(Set2))
  
  cluster.match = list(Set1= rep('',n1), 
                       Set2= rep('',n2))
  names(cluster.match[[1]]) <- names(Set1) 
  names(cluster.match[[2]]) <- names(Set2) 

 
  y1 = rep(0, n1)
  for(i in 1:n1){
    tmp <- lapply(Set2, function(x) jaccard.sim (Set1[[i]], x)) %>% unlist() 
    cluster.match[[1]][i] = which.max(tmp) %>% names()
    y1[i] <- tmp[which.max(tmp)] 
  }
  names(y1) <- names(Set1)
  # merge the score for the same CT state
  if(merge.cluser & n1 > length(cluster1)) {
    y1.bk <- y1
    y1 <- NULL
    for(k in 1:length(cluster1)){
      y1 <-c(y1, max(y1.bk[cluster1[k]]))
    }
    n1 <- length(cluster1)
  }
  Recover <- sum(y1)/n1
 
  y2 = rep(0, n2)
  for(i in 1:n2){
    tmp <- lapply(Set1, function(x) jaccard.sim(x, Set2[[i]])) %>% unlist() 
    cluster.match[[2]][i] = which.max(tmp) %>% names()
    y2[i] <- tmp[which.max(tmp)]  
  }
  names(y2) <- names(Set2)
  # merge the score for the same CT state
  if(merge.cluser & n2 > length(cluster2)) {
    y2.bk <- y2
    y2 <- NULL
    for(k in 1:length(cluster2)){
      y2 <-c(y2, max(y2.bk[cluster2[k]]))
    }
    n2 <- length(cluster2)
  }
  Relevance <- sum(y2)/n2
  
  if(weight){
    w1 <- w.score4cluster.match( cluster.match$Set1 )
    w2 <- w.score4cluster.match( cluster.match$Set2 )
    Recover = Recover*w1
    Relevance = Relevance*w2
  }
  
  
  F1 = 2/((1/Recover) + (1/Relevance))
 

  return(list(F1=F1, cluster.match =cluster.match))
  
}


# transform = function(x){ 
#   y = ifelse(x, 1, -1)
#   return(y)
# }
w.score4cluster.match <- function(Set){
  y = Set==names(Set)
  #z <- transform(y)
  z = length(which(y))/length(Set)
  return(z)  
}
  
Normalize.F1 <- function(F1){
  if(all(F1==0) | all(F1==1)) F1.norm = F1 else { 
    F1.norm = pnorm(scale(F1, center = TRUE, scale = TRUE))
  }
  return(F1.norm)
}

AggregatedMean.F1 <- function(F1.norm){
  S = mean(F1.norm)
  return(S)
}

myplotIc <- function(filename,BioTIP_scores, CTS.candidate,SimResults_g ,
                     width=10, height=6, nn=NULL){
  pdf(file=filename, width=width, height= height)
  n = length(BioTIP_scores)
  # x.row= ifelse(n>=4, 2, 1)
  # if(n>8) warning('Only up to 8 CTSs are plotted')
  if(is.null(nn)) nn<- ifelse(n>8, 8 ,n)
  par(mfrow=c(2,nn)) # 8 plots per page, SimResults_g
  for(i in 1:nn){
    x = length(CTS.candidate[[i]])
    plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
                       ylim=c(0, max(unlist(BioTIP_scores))),
                       las = 2, ylab="Ic.shrink",
                       main= paste("Cluster",names(BioTIP_scores)[i],"\n",x," genes"),
                       fun="boxplot",  #fun="matplot", 
                       which2point=names(BioTIP_scores)[i])
    TEXT = ifelse(local.IC.p, "p.Local=", "p.Global=")
    text(1, 0.09,  paste(TEXT,p.IC[i]), cex=1.5 )
  }
  for(i in 1:nn){
    x = length(CTS.candidate[[i]])
    plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
                       main = paste("Delta Ic.shrink",x,"genes"), 
                       ylab=NULL,
                       xlim=range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]],SimResults_g[[i]])))
  }
  
  if(n>8) {
    for(i in (nn+1):n){
      x = length(CTS.candidate[[i]])
      plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
                         ylim=c(0, max(unlist(BioTIP_scores))),
                         las = 2, ylab="Ic.shrink",
                         main= paste("Cluster",names(BioTIP_scores)[i],"\n",x," genes"),
                         fun="boxplot",  #fun="matplot", 
                         which2point=names(BioTIP_scores)[i])
      TEXT = ifelse(local.IC.p, "p.Local=", "p.Global=")
      text(1, 0.09,  paste(TEXT,p.IC[i]), cex=1.5 )
    }
    for(i in (nn+1):n){
      x = length(CTS.candidate[[i]])
      plot_SS_Simulation(BioTIP_scores[[i]], SimResults_g[[i]], 
                         main = paste("Delta Ic.shrink",x,"genes"), 
                         ylab=NULL,
                         xlim=range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]],SimResults_g[[i]])))
    }
    
  }
  dev.off() 
}



get.downsized.F1.score <- function(CT.ref, getTopMCI.gene.minsize, MCIbottom, i=1)
{
  F1 <- F1.ctl <- N.module <- list()
  for(m in 1:length(getTopMCI.gene.minsize)){
    load(file=paste0('stability/BioTIP_res_10runs_MCIbottom',MCIbottom[i],
                     '_Modulesize',getTopMCI.gene.minsize[m],'.RData'))
    F1.CTS <- F1.random.CTS <- n.module <- array(dim=n.scability)
    for(j in 1:n.scability){
      # evaluate the F1 score for two CTS identificitons
      set2 <- res[[j]]$CTS.candidate[which(res[[j]]$significant)] 
      F1.CTS[j] <- F_score.CTS(CT.ref, set2, weight=TRUE)$F1
      # negative control
      n <- lengths(set2)
      random.set2 <- lapply(n, function(x) sample(rownames(sce), x))
      F1.random.CTS[j] <- F_score.CTS(CT.ref, random.set2, weight=TRUE)$F1
      n.module[j] <- n[names(CT.ref)]
    } 
    F1[[m]] =  F1.CTS
    F1.ctl[[m]] =  F1.random.CTS
    N.module[[m]] = n.module
  }
  names(F1) <- names(F1.ctl) <- names(N.module) <- getTopMCI.gene.minsize
  return(list(F1=F1, F1.ctl=F1.ctl, N.module=N.module))
}

get.multiple.F1.score <- function(fid=NULL, match.col.name=NULL,
                                  #getTopMCI.gene.minsize=30, MCIbottom=2,
                                  CTS.reference = NULL,
                                  QuanTC.genes.list.id = 5,
                                  running.methods=c("./C_Leiden_0.4" , "./C_Leiden_0.8" , "./C_Leiden_1.2"), 
                                  methods=c( "Leiden_0.4" ,"Leiden_0.8", "Leiden_1.2" ),
                                  sce,
                                  best.matched.clusterID)
# multiple CTS.reference read from fid
# otherwise useing the CTS.reference (a vector of gene names) for one CT  state one by one
{
  F1.scores <- F1.ctl <- list()
  if(is.null(match.col.name) & is.null(CTS.reference)) stop('pls give a value for eigher match.col.name or CTS.reference')
  # one CTS.reference
  if(!is.null(CTS.reference)) {
    if(!is.null(fid)) warning('using the CTS.reference and ignoring fid')
    fid = 'An_identification'
  }
  # multiple CTS.reference read from fid, need further debug yet !!!!!!!!!!!!!!!!!
  for(m in 1:length(fid)){
    if(is.null(CTS.reference)) {
      CTS.reference <- read.table(fid[m])
      CTS.reference <- CTS.reference[,1]
    }
    
    F1.CTS <- F1.random.CTS <- array()
    for(i in 1:length(running.methods)){
      set1 <- list()
      set1[[1]] = CTS.reference  
      names(set1) <- match.col.name[m] 
      x <- methods[i]
      if(grepl('QuanTC', methods[i] )) {
        load(file=paste0(running.methods[i],'/CTS.RData'))
        # only the transition genes for the early HEP were recorded for comparision
        names(CTS)
        #[1] "C1-C4" "C4-C2" "C6-C5" "C5-C8" "C8-C3" "C7"
        set2 <- CTS[QuanTC.genes.list.id] # the transition involving softer-clustering C3 !!!!!!
        if(length(QuanTC.genes.list.id)>1) set2 <- list(unlist(set2))
        names(set2) <- names(set1)
        F1.CTS[i] <- F_score.CTS(set1, set2, weight=TRUE, merge.cluser=TRUE)$F1
        # negative control
        n <- lengths(set2)
        random.set2 <- lapply(n, function(x) sample(rownames(sce), x))
        F1.random.CTS[i] <- F_score.CTS(set1, random.set2, weight=TRUE)$F1
      } else {
        load(file=paste0(running.methods[i],'/BioTIP.res.RData'))
        if(length(res)>1){
          names(res)
          #[1] "CTS.candidate" "CTS.score"     "Ic.shrink"     "significant" 
          set2 <- res$CTS.candidate[which(res$significant)]
          if(length(set2) >0) {
            matched.name <- best.matched.clusterID[paste0('C_',methods[i]), match.col.name[m] ]  #!!!!
            names(matched.name) <- names(set1)
            names(set1) <- matched.name[names(set1)]
            F1.CTS[i] <- F_score.CTS(set1, set2, weight=TRUE)$F1
            # negative control
            n <- lengths(set2)
            random.set2 <- lapply(n, function(x) sample(rownames(sce), x))
            F1.random.CTS[i] <- F_score.CTS(set1, random.set2, weight=TRUE)$F1
            
          } else  F1.CTS[i] = F1.random.CTS[i] = 0
        } else  F1.CTS[i] = F1.random.CTS[i] = 0
      }
    }  
    names(F1.CTS) <- names(F1.random.CTS) <-methods
    F1.scores[[m]] <- F1.CTS
    F1.ctl[[m]] <- F1.random.CTS
  }
  names(F1.scores) <- names(F1.ctl) <- match.col.name
  return(list(F1.scores=F1.scores, 
              F1.ctl=F1.ctl))
}


find.proxy.GS <- function(n.select, select,best.matched.clusterID, match.col.name='MEP'){
  if(class(best.matched.clusterID)!='data.frame') stop('pls input best.matched.clusterID as dta.frame')
  if(! match.col.name %in% colnames(best.matched.clusterID)) stop ('best.matched.clusterID has no match.col.name')
  GS.CTS = NULL
  for(i in 1:length(select)){
    if(select[i] == 'QuanTC_run') {
      load(file=paste0(select[i],'/CTS.RData'))
      set2 <- CTS[[2]]  # CTS[[1]] for E18.5
      GS.CTS = c(GS.CTS, unlist(set2[x]))
      rm(CTS)
    } else  {
      load(file=paste0('C_',select[i],'/BioTIP.res.RData'))
      set2 <- res$CTS.candidate[which(res$significant)]
      best.match <- best.matched.clusterID[paste0('C_',select[i]),match.col.name]  #!!!!!!!
      x <- grep(best.match, names(set2))
      GS.CTS = c(GS.CTS, unlist(set2[x]))
    }
  }
  x <- table(GS.CTS)
  GS.CTS <- names(x)[which(x>= n.select)]
  return(GS.CTS)
}
