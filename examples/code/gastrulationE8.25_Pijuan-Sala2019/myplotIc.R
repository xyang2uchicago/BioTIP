
myplotIc <- function(filename,BioTIP_scores, CTS.candidate,SimResults_g ,
                     width=10, height=6, local.IC.p = TRUE, nn=NULL){
  require(BioTIP)
  
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

