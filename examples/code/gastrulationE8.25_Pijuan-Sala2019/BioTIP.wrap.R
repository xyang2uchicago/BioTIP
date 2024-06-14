################################################################
#  A wrapped function to run BioTIP 
#################################################################
# Last update 1/14/2022
# Author: Holly Yang (xyang2.at.uchicago.edu)
#' @title A master function to run BioTIP analysis 
#'
#' @description
#' There are two purposes of the \code{BioTIP}() function.  
#' One is to make the tool easy to use for every R user. Another purpose is to enable iteratively calls to estimate the stability of predictions.
#' We recommend expert- or interested users to perform each step, with their own tweaking, as needed.
#'
#' @param  sce a \code(SingleCellExperiment} Object or \code{cell_data_set} Object
#' @param  samplesL A list of vectors, whose length is the number of states. Each vector gives the sample names in a state. 
#' Note that the vector s (sample names) has to be among the column names of the R object 'df'. This is also inherited from \link[graphics]{matplot}.
#' @param  subDir: A string setting the sub path to save output and intermediate results (when verbose=TRUE) and plots (when plot=TRUE).
#'
#' @param  smallest.population.size: An integer, setting the threshold to exclude the smallest cell cluster with less cells  
#' @param  globa.HVG.select A logical value indicate weather to perform HVG selection across the whole dataset. 
#  By defaule, this value is FALSE when the input sce object is already an object for the HVG genes, 
#' @param  dec.pois A DataFrame is returned where each row corresponds to a gene in x (or in subset.row, if specified). 
#'  This is the output of \code{scran::modelGeneVarByPoisson}.
#' @param  n.getTopHVGs An integer setting the number of global HVG. It only functions when globa.HVG.select=TRUE. 
#'  For BioTIP analysis, this number is suggested to be at least two fold larger than that of Differentially expressed genes.
#'  We recommend to use a number larger than 4000.
#'
#' @param  getTopMCI.n.states An integer determines how many modules are evaluated (the top n MCI scores would be evaluated
#' for downstream analysis). It is inherited from the parameter \code{n} of the function \code{\link{getTopMCI}}. 
#' Default is the number of cell clusters tested.
#' @param  getTopMCI.gene.minsize  A numerical value of the minimum number of transcripts in a cluster. A cutoff that determines
#' the smallest cluster in transcripts numbers for downstream analysis. It is inherited from the parameter \code{min} of the function \code{getTopMCI}.
#' @param  getTopMCI.gene.maxsize  A numerical value of the maximum number of transcripts in a cluster. A cutoff that determines
#' the largest cluster in transcripts numbers for downstream analysis. Default is NULL. 
#' When there are cell poulations over 1k cells, recommend to set this parameter (e.g., 100) to reduce false positive.
#' @param  MCIbottom: a numeric value setting the threshold for MCI score (i.e., the CTS score in the publication), 
#'            below which the identification will be removed for the downstream investigation.
#'
#' @param  local.HVG.optimize A logical value indicating wheather to run function code{\link{optimize.sd_selection}}.
#' Defaul is TRUE. However, for single-cell RT-PCR data or dataset with small number of tested geens, we recommend to set this parameter to FALSE, i.e., 
#' to run BioTIP analysis on all global HVG genes.
#' @param  logmat.local.HVG.testres A data matrix of logcounts for the pre-selected local HVGs, 
#' to speed up the process of large datasets by skipping the step of selecting local HVGs.  
#' Default is NULL to select local HVGs iteratively by \code{localHVG.runs} runs, if \code{local.HVG.optimize} is TRUE. 
#' @param  localHVG.preselect.cut A positive numeric between 0 and 1, setting how many propotion of global HVG to be selected per cell cluster.
#' Default is 0.01. See the parameter \code{cutoff} of the functions code{\link{optimize.sd_selection} and/or code{\link{sd_selection}.
#' @param  localHVG.runs An integer indicating the number of times to run this optimization, default 1000.
#'
#' @param  getNetwork.cut.fdr A numeric cutoff value for a Pearson Correlation Coefficient (PCC) analysis. 
#' Default is 0.05. Transcripts are linked into a network if their correlations meet this PCC-significance criterion
#' See the parameter \code{fdr} of the function \code{\link{getNetwork}} for details. It only functions when code{local.HVG.optimize} is TRUE.
#' @param  n.getMaxMCImember an integer deciding to report n modules with leading MCI scores per state (cluster).
#'   Default is 2, i.e., two gene modules per cluster will be evaluated for significance.  
#'   It inherits the parameter \code{n} of the function \code{\link{getMaxMCImember}}.
#' @param  getMCI.adjust.size A Boolean value indicating if MCI score should be adjusted by module size (the number of transcripts in the module) or not. Default FALSE. 
#' This parameter is not recommended for fun=BioTIP. This parameter is called by the function \code{getMCI}.
#' 
#' @param  n.permutation: An integer setting the number of permutation runs to estimate the significance of CTS identifications
#' @param  empirical.MCI.p.cut  A value setting the threshold of significance when estimating the significance of CTS identifications using the DNB model 
#' 
#' @param  empirical.IC.p.cut  A value setting the threshold of significance when estimating the significance of CTS identifications using the Ic model
#' @local.IC.p  A logical value. TRUE indicates that the singificance comes by comparing the observed Ic to simulated Ic of the same cluster;
#' FLASE will comapre the observed Ic to the simulated Ic at all tested clusters. 
#' @IC.rank  An integer, indicating that the singnificance comes also by observing the rank of the IC score at the desired cell cluster
#' among the IC scores for all tested cell clusters in the dataset. Default is 1 to request the highest rank. For complex datasets, this parameter can be adapted to 
#' a larger number. 
#' 
#' @param  M is precalculated gene-gene correlation matrix, will be reused in each downstream simulation analysis
#   M is also reusable when rerun on the same set of genes and samples. 
#' @param  permutation.method An option of the way to estimate the significance of CTS identifications. 
#' For simple trajectory topology or simple dataset focusing on one branch, we recommend 'both', i.e., shuffling both genes and cells.
#' For complex or large-ranged trajectory with many branches, we recommend 'gene', i.e., only shuffling genes
#' @param  verbose A logical value. When TRUE, save intermediate results into the code{subDir}.
#' @param  plot A logical value. When TRUE, generate plots for the intermediate results into the \code{subDir}.
#'
#'
#' @return A list of four objects when there is an identification. This contains the fields:
#' @param   CTS.candidate is a named list of identified co-expressed gene models (i.e., CTS candidate genes) where the name indicates the predicted cell cluster IDs.
#' @param   CTS.score is a named vector of numbers, where the name indicates the predicted cell cluster IDs.  
#' @param   Ic.shrink is a named list of Ic.shrink scores with the name of this list indicating the predicted cell cluster ID based on a CTS.candidate.
#' Each element of the list is a vector of Ic.shrink scores calculated across all input cell clusters for the corresponding CTS.candidate. 
#' @param   significant A vector of named logical values indicating whether a CTS.candidate at one cluster (the name) is significant in both models of DNB and Ic.shrink.
#' 
#' @export

BioTIP.wrap <- function(sce, samplesL, subDir = 'newrun', 
                   smallest.population.size=20,  
                   globa.HVG.select=FALSE,
                   dec.pois=NULL, 
                   n.getTopHVGs=round(nrow(sce)*0.4),  # used for global HGV selection, dependening on the number of genes
                   getTopMCI.n.states=ncol(sce), 
                   getTopMCI.gene.minsize=30, 
                   getTopMCI.gene.maxsize=NULL, 
                   MCIbottom=2,
                   local.HVG.optimize=TRUE, 
                   localHVG.preselect.cut=0.1, 
                   localHVG.runs=100, 
                   logmat.local.HVG.testres = NULL,
                   getNetwork.cut.fdr=0.05, 
                   n.getMaxMCImember = 2,   
                   getMCI.adjust.size =FALSE,
                   n.permutation = 100, 
                   empirical.MCI.p.cut=0.01, 
                   empirical.IC.p.cut = 0.05, 
                   local.IC.p = TRUE, 
                   IC.rank=1,
                   M=NULL,
                   permutation.method=c('gene','both','sample'),
                   verbose=FALSE, plot=TRUE )
{
  require(SingleCellExperiment)
  if(!'logcounts' %in% assayNames(sce) & !'counts' %in% assayNames(sce)) stop("no 'counts' nor logcounts' in names(assays(sce))")
  if(! permutation.method %in% c('gene','both','sample')) {
    permutation.method='gene'
    warning('No permutation.method found for Ic signicance, gene permutation is performed')
  }
  require(ggplot2)
  
  mainDir = getwd()
  if (!file.exists(subDir)){
    dir.create(file.path(mainDir, subDir))
   # setwd(file.path(mainDir, subDir))
  }
  outputpath = paste0(subDir,'/BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,"_")
  
  (tmp = lengths(samplesL))
  if(any(tmp < smallest.population.size))  samplesL <- samplesL[-which(tmp < smallest.population.size)]
  sce = sce[,unlist(samplesL)]
  
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
    
     
  ######## 1) select global HVG ######################
  if(globa.HVG.select){
    if(is.null(dec.pois)){
     ## set.seed(0010101)
      dec.pois <- scran::modelGeneVarByPoisson(sce,  block=sce$batch)
      } else dec.pois = dec.pois[rownames(sce),]
       
    libsf.var <- getTopHVGs(dec.pois, n=n.getTopHVGs)
    length(libsf.var)
    # [1] 4000
    table(libsf.var %in% rownames(sce))
    # FALSE  TRUE
    #  927  3073
    libsf.var <- intersect(libsf.var, rownames(sce))
    length(libsf.var) # 3073
    dat <- sce[libsf.var,]
  } else dat <- sce
  
  #logcounts normalized by library size
  if(grepl("cell_data_set", class(dat))) {
    require(monocle3)
    logmat <- as.matrix(normalized_counts(dat))
  } else if((grepl("SingleCellExperiment", class(dat)))) logmat <- as.matrix(logcounts(dat))
  
  dim(logmat)
  # [1] 3073 7240
  rm(dat)
  ############ 2) identify CTS candidates 
  ### BEGIN DO NOT REPEAT #1 -----------------------------
  # 
  # 2.1)  select local HVG per cluster by this optimizes sd selection
  if(local.HVG.optimize) {
    if(is.null(logmat.local.HVG.testres)){  
   ## set.seed(2020)
    cat('running local HVG optimalization ...')
    # testres is a list of downsampled logmat per cluster, with rows to be locally selected HVGs
    testres <- optimize.sd_selection(logmat, samplesL, B=localHVG.runs, 
                                   cutoff=localHVG.preselect.cut,
                                   times=.75, percent=0.8)
    if(verbose) save(testres, file=paste0(outputpath,"optimized_local.HVG_selection.RData"), compress=TRUE) # !!!!!!!!!!!!!!!!
    } else  {
      testres <- logmat.local.HVG.testres # working on the pre-selected local HVGs
      if(!all(rownames(testres) %in%  rownames(logmat))) {
        testres <- testres[intersect(rownames(testres),rownames(logmat)),]
        warning('Some of genes in the given logmat.local.HVG.testres are outside the golbal HVG!')
      }
    }
  }  else {
    testres <- sd_selection(logmat,  samplesL,  cutoff=localHVG.preselect.cut) 
  }
  # load(file=paste0(outputpath,"optimized_local.HVG_selection.RData"))
  names(testres)
  cat('HVG selection is done', '\n')
  
  # Network Partition
  igraphL <- getNetwork(testres, fdr = getNetwork.cut.fdr)
  # C1:257 nodes
  # C2:264 nodes
  # C3:286 nodes
  # C4:295 nodes
  # TC:211 nodes
  
  cluster <- getCluster_methods(igraphL)
  names(cluster)
  # [1] "C1" "C2" "C3" "C4" "TC"
  cat('Network partition is done', '\n')
  
  # 2.2) Identifying CTS candidate, the MCI score uisng the DNB model  #########
  membersL <- getMCI(cluster,testres, adjust.size = getMCI.adjust.size, fun='BioTIP')
  names(membersL)
  #[1] "members" "MCI"     "sd"      "PCC"     "PCCo"
  # save(membersL, file=paste0(outputpath,"membersL.RData")) #!!!!!!!!!!!!!!
  cat('CTS score calculation is done', '\n') 
  
  # get the statistics using the MCI system
  # getTopMCI.n.states decides how many states to check
  topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], 
                     membersL[["MCI"]], min=getTopMCI.gene.minsize, 
                     n=getTopMCI.n.states)
  topMCI
  #   C2       C1       TC 
  # 4.209349 3.518579 3.413213 
  
  if(class(topMCI)=="numeric" & length(topMCI)>0){
    if(plot){
      w <- ifelse(length(membersL$MCI) > 10, 80, 40)
      h <- ifelse(length(membersL$MCI) > 10, 10, 5)
      pdf(file=paste0(outputpath,"MCIBar_", localHVG.preselect.cut,
                      "_fdr",getNetwork.cut.fdr,"_minsize",getTopMCI.gene.minsize,".pdf"),
          width=w, height=h)
      plotBar_MCI(membersL, ylim=c(0, ceiling(max(topMCI, na.rm=TRUE))), minsize = getTopMCI.gene.minsize)
      if(!is.null(MCIbottom)) abline(h=MCIbottom, lty=2, col='red')
      plotBar_MCI(membersL, ylim=c(0, ceiling(max(topMCI, na.rm=TRUE))*2))
      if(!is.null(MCIbottom)) abline(h=MCIbottom, lty=2, col='red')
      dev.off()
    }
    # #the numbers in the parenthesis: they are total number of clusters, no control of the cell (sample) size per  cluster
    
    
    x <- which(topMCI >= MCIbottom)
    if(length(x)>0) {
      topMCI = topMCI[x] 
      if(getTopMCI.n.states > length(x)) warning(paste('less number of states have the highest CTS score larger than',MCIbottom))
      getTopMCI.n.states = length(x)
      
      # get the CTS ID within each cluster and CTS.candidate genes for the first leading MCI scores per state
      CTS.candidate.ms <- getMaxMCImember(membersL[["members"]],
                                          membersL[["MCI"]],min =getTopMCI.gene.minsize, 
                                          n=n.getMaxMCImember)
      names(CTS.candidate.ms)
      CTS.candidate = getCTS(topMCI, CTS.candidate.ms[["members"]][names(topMCI)])
      # Length: 54
      # Length: 26
      # Length: 51
      
      if(n.getMaxMCImember>1){
        # get the state ID and CTS.candidate genes for the 2nd leading MCI scores per state
        ## tmp extracts the number of bars within each named state
        tmp <- unlist(lapply(CTS.candidate.ms[['idx']][names(topMCI)], length))
        ## here returns all the groups with only 2 or more bars
        (whoistop2nd <- names(tmp[tmp>=2]))
        if(length(whoistop2nd)>0) {
          nextMCI = getNextMaxStats(membersL[['MCI']], idL=CTS.candidate.ms[['idx']], whoistop2nd, which.next = 2) 
          nextMCI = nextMCI[order(nextMCI, decreasing=TRUE)]
          x <- which(nextMCI >= MCIbottom)
          if(length(x)>0) { 
            nextMCI = nextMCI[x]
            whoistop2nd = names(nextMCI)
            ## add the gene members of the 2nd toppest biomodue in the states with exactly 2 bars
            CTS.candidate = append(CTS.candidate, CTS.candidate.ms[["2topest.members"]][whoistop2nd])
            topMCI = append(topMCI, nextMCI)
            ## here returns all the groups with exactly 3 bars
            if(n.getMaxMCImember>2){
              whoistop3rd = names(tmp[tmp>=3])
              if(length(whoistop3rd)>0)  {
                nextMCI = getNextMaxStats(membersL[['MCI']], idL=CTS.candidate.ms[['idx']], whoistop3rd, which.next = 3) 
                nextMCI = nextMCI[order(nextMCI, decreasing=TRUE)]
                x <- which(nextMCI >= MCIbottom)
                if(length(x)>0) { 
                  CTS.candidate = append(CTS.candidate, CTS.candidate.ms[["3topest.members"]][whoistop3rd])
                  topMCI = append(topMCI, nextMCI) 
                  if(n.getMaxMCImember>3) warning('Please manually modify the BioTIP() accordingly to extract 4th and more CTS candididate per cell cluster')
                }
              }  
            }
          }   
        } # ending  if(length(whoistop2nd)>0)    
      } else warning(paste('All 2nd highest CTS scores are lower than',MCIbottom))
      
      if(verbose) save(CTS.candidate, topMCI, file=paste0(outputpath,"CTS.candidate.RData")) #!!!!!!!!
      rm(sce)
      ### END DO NOT REPEAT, #1 -----------------------------
      
      
      # load(file=paste0('BioTIP_top',localHVG.preselect.cut,'FDR',getNetwork.cut.fdr,'_','CTS.candidate.RData'))
      # names(CTS.candidate)
      # # [1] "C2" "C1" "TC"
      
      
      ############# estimate significance for CTS cendidate using the DNB model 
      # ## BEGIN DO NOT REPEAT, #2----------------------------- 
      # also refer to 'BioTIP_E8.25_mesoderm_cluster.R' 
      
      dim(logmat)
      #[1]  3073 1362
      # M is precalculated gene-gene correlation matrix, will be reused in the downstream simulation analysis
      if(is.null(M)) {
        
        cat('calcualting M .... ')
        M <- cor.shrink(logmat, Y = NULL, MARGIN = 1, shrink = TRUE)
        save(M, file=paste0(outputpath,"CTS_ShrinkM.RData"), compress=TRUE) # save the file as its calculation runs a while
        
      }
      dim(M)  #3073  3073
      M = M[rownames(logmat), rownames(logmat)]
      
      # load(file=paste0(outputpath,"CTS_ShrinkM.RData"))
      
      
      simuMCI = list()
      ## set.seed(2020)
      for (i in 1:length(CTS.candidate)){
        n <- length(CTS.candidate[[i]])
        simuMCI[[i]] <- simulationMCI(n, samplesL, logmat,  adjust.size=getMCI.adjust.size, B=n.permutation, fun="BioTIP", M=M)
      }
      names(simuMCI) = names(CTS.candidate)
      
      if(verbose) save(simuMCI, file=paste0(outputpath,"SimuMCI_",n.permutation,"_", localHVG.preselect.cut, "_fdr",
                                            getNetwork.cut.fdr,"_minsize",getTopMCI.gene.minsize,".RData"))
      # 
      
      n = length(CTS.candidate)
      if(plot & n>0 & class(topMCI)=="numeric" & length(topMCI)>0){
        w <- ifelse(n>3, 4*n, 7); h <- ifelse(n>3, n, 7) # updated 1/28/2022
        pdf(file=paste0(outputpath,"barplot_MCI_Sim_RandomGene.pdf"), width=w, height=h)
        par(mfrow=c(1,n))  
        for (i in 1:n){
          plot_MCI_Simulation(topMCI[i], simuMCI[[i]], las=2, ylim=c(0, max(c(topMCI, simuMCI[[i]]), na.rm=TRUE)),
                              main=paste("Cluster", names(CTS.candidate)[i], ";", length(CTS.candidate[[i]]), "genes",
                                         "\n","vs. ",n.permutation, "times of gene-permutation"),
                              which2point=names(CTS.candidate)[i])
        }
        dev.off()
      }
      
      ## END DO NOT REPEAT #2 -----------------------------
      
      ## dropoff non-significant CTS by the above MCI simulation ---------------
      
      n <- length(CTS.candidate)
      p <- array(dim=n)
      for(i in 1:n){
        ## scan the global MCI score in all clusters:
        p[i] = length(which(simuMCI[[i]] >= topMCI[i]))/n.permutation/length(samplesL) # updated 2/4/2022
      }
      dropoff <- which(p>= empirical.MCI.p.cut)  
      if(any(dropoff)){
        CTS.candidate <- CTS.candidate[-dropoff]  # updated 1/29/2022
        #  CTS.Symbol <- CTS.Symbol[-dropoff]
      } 
      
      # save(CTS, file= paste0(outputpath, 'CTS.candidate.RData'))
      # rm(n, CTS.candidate)
      # load(file= paste0(outputpath, 'CTS.candidae.RData'))
      cat('CTS.dandidate is done', '\n')
      
      ########  Finding Tipping Point and verify using IC score !!!!  #################
      ######  BioTIP score, evaluate by Ic_shrink scores by shuffling genes, samples, or both ##############
      if(length(CTS.candidate)>0){
        # C= 1000
        BioTIP_scores <- SimResults_g <- SimResults_s <-  SimResults_b <- list()
        ## set.seed(101010)
        
        for(i in 1:length(CTS.candidate)){ # begin loop i
          CTS <- CTS.candidate[[i]]
          n <- length(CTS)
          BioTIP_scores[[i]] <- getIc(logmat, samplesL, CTS, fun="BioTIP", 
                                      shrink=TRUE, PCC_sample.target = 'none' )
          if(permutation.method=='both'){ # flag1
            ######  shuffling both genes and samples ##############
            SimResults_b[[i]] <- matrix(nrow=length(samplesL), ncol=n.permutation)
            rownames(SimResults_b[[i]]) = names(samplesL)
            for(j in 1:length(samplesL)) {
              ns <- length(samplesL[[j]])  # of each state
              CTS.sim <- sample(rownames(logmat), n) # !!!!
              SimResults_b[[i]][j,] <- simulation_Ic_sample(logmat, ns, #Ic=BioTIP_score[i],
                                                            genes= CTS.sim, B=n.permutation,
                                                            fun="BioTIP", shrink=TRUE, PCC_sample.target = 'none')
            } 
          } else { # flag1.1
            if(permutation.method=='sample'){ # flag2
              #### shuffling cluser IDs for cells
              SimResults_s[[i]] <- matrix(nrow=length(samplesL), ncol=n.permutation)
              rownames(SimResults_s[[i]]) = names(samplesL)
              for(j in 1:length(samplesL)) {
                ns <- length(samplesL[[j]])  # of each state
                SimResults_s[[i]][j,] <- simulation_Ic_sample(logmat, ns, #Ic=BioTIP_score[i],
                                                              genes=CTS, B=n.permutation,
                                                              fun="BioTIP", shrink=TRUE, PCC_sample.target = 'none')
              }
            } else { # flag2.1
              #### shuffling genes
              SimResults_g[[i]]  <- simulation_Ic(n, samplesL, logmat, B=n.permutation,
                                                  fun="BioTIP", shrink=TRUE, PCC_sample.target = 'none')
              
            } # flag2.2
          } # flag1.2
          
        } # end loop i
        
        # names(SimResults_b) <- names(CTS.candidate)
        
        names(BioTIP_scores) <- names(CTS.candidate)
        if(length(SimResults_g)>0) names(SimResults_g) <- names(BioTIP_scores) 
        if(length(SimResults_s)>0) names(SimResults_s) <- names(BioTIP_scores) 
        if(length(SimResults_b)>0) names(SimResults_b) <- names(BioTIP_scores) 
        
        
        if(verbose) if(permutation.method=='both') {
          save( SimResults_b, BioTIP_scores, 
                file=paste0(outputpath,"IC_sim.PermutateBoth.RData"), compress=TRUE)
        } else {
          save( SimResults_g, SimResults_s, BioTIP_scores, 
                file=paste0(outputpath,"IC_sim.Permutation.RData"), compress=TRUE)
        }
        
        # 
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
            interesting =  names(BioTIP_scores[i]) ## uodated 1/30/2022
            #first p value calculated for exactly at tipping point
            p = length(which(SimResults_4_p.IC[[i]][interesting,] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
            p = p/n.permutation  # ncol(SimResults_4_p.IC[[i]])
            #second p value calculated across all statuses
            p2 = length(which(SimResults_4_p.IC[[i]] >= BioTIP_scores[[i]][names(BioTIP_scores)[i]]))
            p2 = p2/n.permutation  #ncol(SimResults_4_p.IC[[i]])
            p.IC[i] = p2/nrow(SimResults_4_p.IC[[i]])
          }
        }   
        cat('BioTIP is done', '\n')
        
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
          if (IC.rank==1) x2[i] = ifelse(names(which.max(BioTIP_scores[[i]])) == names(BioTIP_scores)[i],
                                         TRUE, FALSE) else{
                                           n <- length(BioTIP_scores[[i]])
                                           which.high <- rank(BioTIP_scores[[i]])[n:(n-IC.rank+1)] 
                                           x2[i][which.high] <- TRUE
                                         }
        }
        significant = x & x2 
        names(significant) <- names(CTS.candidate)
        
      } else { 
        BioTIP_scores = NULL
        significant = rep(FALSE, length(CTS.candidate))
        names(significant) <- names(CTS.candidate)
        cat('All CTS.candidate are insignificant in DMB model')
      }
      
      if(!is.null(getTopMCI.gene.maxsiz)) {
        x <- lengths(CTS.candidate)
        dropoff <- which(x> getTopMCI.gene.maxsize)
        if(length(dropoff)>0){
          CTS.candidate = CTS.candidate[-dropoff]
          topMCI = topMCI[-dropoff]
          BioTIP_scores = BioTIP_scores[-dropoff]
          significant = significant[-dropoff]
        }
      }
      
      return(list(CTS.candidate=CTS.candidate, 
                  CTS.score=topMCI, 
                  Ic.shrink=BioTIP_scores, 
                  significant=significant)) 
      
    } else {
      warning(paste('No CTS has a score higher than', MCIbottom, '!'))
      return(CTS.candidate=NULL)
    }
  } else{
    warning(paste('No module has', getTopMCI.gene.minsize, ' gene members!'))
     return(CTS.candidate=NULL)
  }
  
 # setwd(file.path(mainDir))
  
}

