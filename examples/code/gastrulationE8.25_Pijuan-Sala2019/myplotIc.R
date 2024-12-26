
myplotIc <- function(filename,BioTIP_scores, CTS.candidate, SimResults_g , p.IC,
                     width=10, height=6, local.IC.p = TRUE, nn=NULL){
  require(BioTIP)

  if(is.null(names(CTS.candidate))) stop('pls give a named list of "CTS.candidate" ')
   ## ensure the same order                       
  if(all(names(BioTIP_scores) %in% names(CTS.candidate))) names(BioTIP_scores) = names(CTS.candidate) else print('names(BioTIP_scores) != names(CTS.candidate)')
  if(all(names(SimResults_g) %in% names(CTS.candidate))) names(SimResults_g) = names(CTS.candidate) else print('names(SimResults_g) != names(CTS.candidate)')
  if(all(names(p.IC) %in% names(CTS.candidate))) names(p.IC) = names(CTS.candidate) else print('names(p.IC) != names(CTS.candidate)')
    
  ## reformate the names of MCI to be within the original names of cell clsuters for visualization ; added by xyang2
  if(any(grepl('.', names(BioTIP_scores), fixed=TRUE))){
    names(BioTIP_scores) = lapply(names(BioTIP_scores), function(x) unlist(strsplit(x, split='.', fixed=T))[1]) %>% unlist
	  names(CTS.candidate) = lapply(names(CTS.candidate), function(x) unlist(strsplit(x, split='.', fixed=T))[1]) %>% unlist
	  names(SimResults_g) = lapply(names(SimResults_g), function(x) unlist(strsplit(x, split='.', fixed=T))[1]) %>% unlist
    names(p.IC) = lapply(names(p.IC), function(x) unlist(strsplit(x, split='.', fixed=T))[1]) %>% unlist
  }
    
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

