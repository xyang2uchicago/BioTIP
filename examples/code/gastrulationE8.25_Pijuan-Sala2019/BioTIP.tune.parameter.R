
BioTIP.tune.parameter <- function(sce, samplesL, subDir = 'newrun', globa.HVG.select=FALSE,
                        dec.pois=NULL, n.getTopHVGs=round(nrow(sce)*0.4),  # used for global HGV selection, dependening on the number of genes
                        getTopMCI.n.states=3, getTopMCI.gene.minsize=30, MCIbottom=2,
                        local.HVG.optimize=TRUE, localHVG.preselect.cut=0.1, localHVG.runs=100,
                        getNetwork.cut.fdr=0.05, 
                        n.getMaxMCImember = 2, 
                        n.permutation = 100, empirical.MCI.p.cut=0.01, 
                        empirical.IC.p.cut = 0.05, local.IC.p = TRUE,
                        M=NULL,
                        permutation.method=c('gene','both','sample'),
                        plot=TRUE)
{
  #require(utils)
  require(ggplot2)
  
  mainDir = getwd()
  if (!file.exists(subDir)){
    dir.create(file.path(mainDir, subDir))
    # setwd(file.path(mainDir, subDir))
  }
  outputpath = paste0(subDir,'/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
  
  load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
  testres
  
  CTS.candidate=NULL
  load(file=paste0(outputpath,"CTS.candidate.RData"))
  if(length(CTS.candidate)>0){
    BioTIP_scores <- SimResults_g <- SimResults_s <-  SimResults_b <- list()
    load(file=paste0(outputpath,"SimuMCI_",n.permutation,"_", localHVG.preselect.cut, "_fdr",
                     getNetwork.cut.fdr,"_minsize",getTopMCI.gene.minsize,".RData"))
    if(permutation.method=='both') {
      load(file=paste0(outputpath,"IC_sim.PermutateBoth.RData"))
    } else {
      load(file=paste0(outputpath,"IC_sim.Permutation.RData"))
    }
    
    n <- length(CTS.candidate)
    p.IC <- rep(1, n)
    if(length(SimResults_g)>0) {
      SimResults_4_p.IC = SimResults_g } else if(length(SimResults_s)>0) {
        SimResults_4_p.IC = SimResults_s 
      } else if(length(SimResults_b)>0) {
        SimResults_4_p.IC = SimResults_b
      }
    if(length(SimResults_4_p.IC)>0) {
      for(i in 1:n){
        interesting = names(BioTIP_scores[i])
        #first p value calculated for exactly at tipping point
        p = length(which(SimResults_4_p.IC[[i]][interesting,] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
        p = p/n.permutation  # ncol(SimResults_4_p.IC[[i]])
        #second p value calculated across all statuses
        p2 = length(which(SimResults_4_p.IC[[i]] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
        p2 = p2/n.permutation  #ncol(SimResults_4_p.IC[[i]])

        if(local.IC.p) p.IC[i] = p/nrow(SimResults_4_p.IC[[i]]) else p.IC[i] = p2/nrow(SimResults_4_p.IC[[i]])
      }
    }   
    
    myplotIc <- function(filename,BioTIP_scores, CTS.candidate,SimResults_g ){
      pdf(file=filename, width=10, height=6)
      n = length(BioTIP_scores)
     # x.row= ifelse(n>=4, 2, 1)
     # if(n>8) warning('Only up to 8 CTSs are plotted')
      nn<- ifelse(n>8, 8 ,n)
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
    
    if(plot){
        if(length(SimResults_g)>0 & length(CTS.candidate)>0) {
          pdf.file.name = paste0(outputpath,"IC_Delta_SimresultGene.pdf")  
          myplotIc(pdf.file.name, BioTIP_scores, CTS.candidate, SimResults_g )
        }
        if(length(SimResults_s)>0) {
          pdf.file.name = paste0(outputpath,"IC_Delta_SimresultSample.pdf") 
          myplotIc(pdf.file.name, BioTIP_scores, CTS.candidate, SimResults_s )
        }
        if(length(SimResults_b)>0) {
          pdf.file.name = paste0(outputpath,"IC_Delta_SimresultBoth.pdf") 
          myplotIc(pdf.file.name, BioTIP_scores, CTS.candidate, SimResults_b )
        }
        
      }
      

    
    ## evaluate by the significance of Ic.shrink after permutation genes  
    n = length(BioTIP_scores)
    # cat('NOTICE!','\n')
    #x <- which(p.IC < empirical.p.cut)
    x <- p.IC < empirical.IC.p.cut  # updated 1/28/2022
    # check if the max Ic is observed at TC
    x2 <- rep(FALSE, n) 
    for(i in 1:n){
      x2[i] = ifelse(names(which.max(BioTIP_scores[[i]])) == names(BioTIP_scores)[i],
                     TRUE, FALSE)
    }
    significant = x & x2 
    names(significant) <- names(CTS.candidate)
    
    return(list(CTS.candidate=CTS.candidate, CTS.score=topMCI, Ic.shrink=BioTIP_scores, significant=significant)) 
    
  } else {
  warning(paste('No CTS has a score higher than', MCIbottom, '!'))
  return(CTS.candidate=NULL)
  }
}
