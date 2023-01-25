## There are four sections in this code:
## section 1) preprocess, gene and cell filtering, and size-factor normalization
## scetion 2) cluster cells (SNNGraph_k5) and build trajectory
## scetion 3) focused on 13 cell clusters, then prepare inputs for QuanTC and BioTIP running, respectively
## section 4) applying cell clustering with different methods

#setwd("F:/projects/scRNA/results/AJ/GSE130146_xy/Results_1k/LibrarySize") 
setwd("scRNAseq_examples/result/EB_Zhao2019")

################################################################################  
## section 1) preprocess, gene and cell filtering, size-factor normalization  ##
################################################################################  

{
  #---------------------------
  # 10x genomics data loading
  # three matrixes were downloaded from 
  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3732839
  # and save in  alocal fold "Data/EB/"
  #---------------------------
  library(DropletUtils)
  fnameeb <- file.path("Data/EB/")
  eb.sce <- read10xCounts(fnameeb, col.names=TRUE)
  
  library(AnnotationHub)
  ens.mm.v97 <- AnnotationHub()[["AH73905"]]
  anno <- select(ens.mm.v97, keys=rownames(eb.sce), 
                 keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))
  rowData(eb.sce) <- anno[match(rownames(eb.sce), anno$GENEID),]
  
  library(AnnotationHub)
  ens.mm.v97 <- AnnotationHub()[["AH73905"]]
  chr.loc <- mapIds(ens.mm.v97, keys=rownames(eb.sce),
                    keytype="GENEID", column="SEQNAME")
  is.mito <- which(chr.loc=="MT")
  #snapshot date: 2020-4-27
  
  # ------------------------------------
  # genes expressed in at least one cells were kept for analysis. 
  # remove noise followed the author's method which were
  # Cells with more than 5% mitochondria reads or fewer than 2000 unique genes detected
  # were filtered out;  
  # we kept 33456 genes and 1731 cells for further analysis. 
  # ------------------------------------
  
  # drop where total values are 0 in df --------------------------
  library(scater)
  df <- perCellQCMetrics(eb.sce, subsets=list(Mito=is.mito))
  df
  
  qc.nexprs <- df$detected < 2000
  qc.mito <- df$subsets_Mito_percent > 5
  discard <- qc.nexprs | qc.mito
  DataFrame(NExprs=sum(qc.nexprs), MitoProp=sum(qc.mito), Total=sum(discard))
  
  filtered <- eb.sce[,!discard]
  dim(filtered)
  #[1] 33456  1731
  
  # save(filtered,file=paste0(mynormalization,"ebfiltered.RData"))
  
  #-----------------------------------
  # normalization
  # The GEO-downloaded (GSE130146) cell expression values were normalized by library size factors.
  #------------------------------------------
  
  library(scater)
  packageVersion('scater') # '1.18.0'
  library(scran)
  packageVersion('scran') # '1.18.0'
  
  set.seed(100)
  filteredsum <- computeSumFactors(filtered, cluster=clust.filtered, min.mean=0.1)
  libsf.sce <- logNormCounts(filteredsum, size_factors=libsf.filtered)
  assayNames(filteredsum)
  
  #save(libsf.sce,file='libsf.sce.RData',compress = TRUE)
  #load('libsf.sce.RData')
  
  libsf.sce
  # class: SingleCellExperiment
  # dim: 33456 1731
  # metadata(1): Samples
  # assays(2): counts logcounts
  # rownames(33456): ENSMUSG00000064842 ENSMUSG00000051951 ...
  # ENSMUSG00000096730 ENSMUSG00000095742
  # rowData names(2): ID Symbol
  # colnames(1731): AAACCTGAGTTTCCTT-1 AAACCTGCATGAAGTA-1 ...
  # TTTGTCATCGGCTTGG-1 TTTGTCATCTGCTTGC-1
  # colData names(3): Sample Barcode sizeFactor
  # reducedDimNames(0):
  #   altExpNames(0):
  #  
}
 
################################################################################  
## section 2) cluster cells and build trajectory,                             ##
################################################################################  
{
  #--------------------------------------------------------------------------
  # dimensionality reduction - run pca first, then TSNE
  #--------------------------------------------------------------------------
  library(scran)
  lib.sce <- modelGeneVar(libsf.sce)
  lib.sce[order(lib.sce$bio, decreasing=TRUE),]
  libsf.var.1000 <- getTopHVGs(lib.sce, n=1000)
  
  library(scater)
  
  set.seed(100) # See below.
  dimred.sce <- runPCA(libsf.sce, subset_row=libsf.var.1000)  # !!!!!!!!!!!!!!!!!
  reducedDimNames(dimred.sce) # 'PCA'
  dim(reducedDim(dimred.sce, "PCA")) # 1731  50
  set.seed(00101001101)
  dimred.sce <- runTSNE(dimred.sce, dimred="PCA") #, ncomponents = 3)
  plotReducedDim(dimred.sce, dimred="TSNE")
  
  #--------------------------------------------------------------------------
  # Clustering cells
  #--------------------------------------------------------------------------
  g.5 <- buildSNNGraph(dimred.sce, k=5, use.dimred = 'PCA')
  clust.5 <- igraph::cluster_walktrap(g.5)$membership
  table(clust.5)
  #clust.5
  #   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17 
  # 102  62 289  98  86  75 151 160 240  38 208  34  19  39  72  47  11  
  
  
  library(scater)
  colLabels(dimred.sce) <- factor(clust.5)
  
  save(dimred.sce, file="GSE130146_robustness/dimred.sce_afterPCA.RData")
  
  
  
  #-----------------------------------------------------------------------------
  # Finding cluster-specific markers
  #-----------------------------------------------------------------------------
  
  library(scran)
  markers.eb <- findMarkers(dimred.sce)
  markers.eb
  
  library(pheatmap)
  
  
  # marker genes for EB from paper except Acta2 and Tnnt2:
  markers <- list()
  markers[['cardiac']] <- c('Isl1','Hand1','Tnnt2')
  markers[['muscle']] <- c('Foxf1','Tbx20','Bmp4','Acta2')
  markers[['hemangiogenic']] <- c('Fli1','Etv2','Tal1','Gata2','Lmo2','Eomes')
  markers[['pluripotent']] <- c('Spp1', 'Sox2','Zfp42','Dppa3','Utf1')
  markers[['early.mesoderm']] <- c('Kdr','Pdgfra','Tbx6','Mesp1','T','Pou5f1')
  markers[['mesoderm']] <- c('Dnmt3b')
  markers[['endoderm']] <- c('Sox17','Fgf5', 'Frzb','Cd24a','Foxa2')
  markers[['cycle']] <- c('Pcna','Top2a', 'Mcm6','Mki67')
  markers[['germ']] <- c('Lefty1','Dnd1')
  all(unlist(markers) %in% rowData(dimred.sce)$Symbol)
  #[1] TRUE
  
  mymarkers <- unlist(markers)
  length(mymarkers) # 36
  
  #save(markers, file="EB_markers.RData")
  
  pdf(file="EB_expression_markers.pdf", width=14)
  for(x in 0:3)
  {
    gridExtra::grid.arrange(
      plotExpression(dimred.sce, features=mymarkers[1+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[1+9*x]),
      plotExpression(dimred.sce, features=mymarkers[2+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[2+9*x]),
      plotExpression(dimred.sce, features=mymarkers[3+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[3+9*x]),
      plotExpression(dimred.sce, features=mymarkers[4+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[4+9*x]),
      plotExpression(dimred.sce, features=mymarkers[5+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[5+9*x]),
      plotExpression(dimred.sce, features=mymarkers[6+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[6+9*x]),
      plotExpression(dimred.sce, features=mymarkers[7+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[7+9*x]),
      plotExpression(dimred.sce, features=mymarkers[8+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[8+9*x]),
      plotExpression(dimred.sce, features=mymarkers[9+9*x], 
                     x="label", colour_by="label",swap_rownames='Symbol',
                     xlab=mymarkers[9+9*x]),
      ncol=3
    )
  }
  
  dev.off()
  
  
  #---------------------------------------------------------------------------
  # trajectory building
  # branching point should be identification
  #-------------------------------------------------------------------------
  
  ## Minimum spanning tree constructed trajectory
  ## The TSCAN package employs a simple yet effective approach to trajectory reconstruction. 
  ## It clusters cells to summarize the data into a smaller set of discrete units, 
  ## computes cluster centroids by averaging the cell coordinates and then forms the minimum spanning tree (MST) 
  ## across centroids. 
  ## The MST is simply an undirected acyclic graph that passes through each centroid exactly once 
  ## and can be thought of as the most parsimonious structure that captures the transitions between clusters.
  
  #load(file="GSE130146_robustness/dimred.sce_afterPCA.RData")
  
  dimred.sce
  # class: SingleCellExperiment 
  # dim: 33456 1731 
  
  library(scater)
  by.cluster <- aggregateAcrossCells(dimred.sce, ids=colData(dimred.sce)$label)
  centroids <- reducedDim(by.cluster, "PCA")
  colData(dimred.sce)
  
  dmat <- dist(centroids)
  dmat <- as.matrix(dmat)
  g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
  mst <- igraph::minimum.spanning.tree(g)
  
  
  pdf(file="trajectory_libsf_1.pdf")
  
  #set.seed(1000)
  plot(mst)
  
  pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
  coords <- reducedDim(by.cluster, "TSNE")
  group <- rep(seq_len(nrow(pairs)), 2)
  stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)
  
  plotTSNE(dimred.sce, colour_by="label", 
           text_by="label", text_size = 8)  
  
  plotTSNE(dimred.sce, colour_by="label", 
           text_by="label", text_size = 8) + 
    geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
  
  
  dev.off()
  
}


################################################################################  
## section 3) focused on 13 cell clusters, then prepare inputs                ##
## for QuanTC and BioTIP running, respectively                                ##
################################################################################  
{
  # After dropping off 4 clusters of endoderm or primordial germ cells, 
  # we kept 1,531 cells along the path from naive pluripotent cells to either hemangiogenic or smooth muscle lineages
  # for further analysis.
  #-------------------------------------------------
  #load(file="GSE130146_robustness/dimred.sce_afterPCA.RData")
  dim(dimred.sce) 
  # dim: 33456 1731 
  
  #dropping endoderm
  dropcols = colData(dimred.sce)$label == 16 |colData(dimred.sce)$label == 2  |
    colData(dimred.sce)$label == 15 |colData(dimred.sce)$label == 13
  dimred.sce = dimred.sce[,!dropcols]
  
  dim(dimred.sce)
  #[1] 33456  1531
  
  
  
  #---------------------------------------------------------------
  # output R object to run BioTIP 
  #---------------------------------------------------------------
  #load(file="GSE130146_robustness/dimred.sce_afterPCA.RData")  
  dim(dimred.sce)
  # dim: 33456 1731 
  
  dropcols = colData(dimred.sce)$label == 16 |colData(dimred.sce)$label == 2  |
    colData(dimred.sce)$label == 15 |colData(dimred.sce)$label == 13
  dimred.sce = dimred.sce[,!dropcols]
  dim(dimred.sce)
  # [1] 33456  1531
  
  colLabels(dimred.sce) <- factor(as.vector(colLabels(dimred.sce)))
  table(colLabels(dimred.sce)) 
  # 1  10  11  12  14  17   3   4   5   6   7   8   9 
  # 102  38 208  34  39  11 289  98  86  75 151 160 240 
  
  # library(scran)
  # lib.sce <- modelGeneVar(logcounts(dimred.sce))
  # lib.sce[order(lib.sce$bio, decreasing=TRUE),]
  # libsf.var <- getTopHVGs(lib.sce, n=4000)
  # save(libsf.var, file='4000_libsf.var_noenderdormPgerm.RData')
  
  # load(file='4000_libsf.var_noenderdormPgerm.RData')
  sce <- dimred.sce[libsf.var,]
  sce
  # class: SingleCellExperiment 
  # dim: 4000 1531 
  # metadata(1): Samples
  # assays(2): counts logcounts
  # rownames(4000): ENSMUSG00000030544 ENSMUSG00000060461 ...
  # ENSMUSG00000070814 ENSMUSG00000031951
  # rowData names(2): ID Symbol
  # colnames(1531): AAACCTGAGTTTCCTT-1 AAACCTGCATGAAGTA-1 ...
  # TTTGTCATCGGCTTGG-1 TTTGTCATCTGCTTGC-1
  # colData names(4): Sample Barcode sizeFactor label
  # reducedDimNames(3): PCA TSNE UMAP
  # altExpNames(0):
  rm(dimred.sce)
  
  #save(sce, file='GSE130146_robustness/sce.GSE130146_noenderdormPgerm.RData')
  save(sce, file='../data/sce.GSE130146_noenderdormPgerm.RData')  # being updated in later sections
  
  
  
  ## cell-cell similarity matrix for 131 cells #---------------------
  # refer to http://127.0.0.1:11637/library/SC3/doc/SC3.html
  # needs to customize ks accordingily per dataset, the larger range the longer running time.
  # In this case, ks=9:19 are tested, 
  # and for the related soft-thresholding clustering (QuanTC method), 
  # we had take average of the Consensus.Cluster-agreeable clustering results of k=4:10 to get cell-cell similarity matrix M
  library(SC3)
  sce.sc3 = sce
  rowData(sce.sc3)$feature_symbol <- rowData(sce.sc3)$Symbol
  # remove features with duplicated names
  if(any(duplicated(rowData(sce.sc3)$feature_symbol))) sce.sc3 <- sce.sc3[-which(duplicated(rowData(sce.sc3)$feature_symbol)),]  # F
  range(assays(sce.sc3)$logcounts)
  # [1]   0.00000 10.59231
  
  ### to run SC3 successfully, transform sparsematrix to matrix  !!!!!!!!!!!!!
  logcounts(sce.sc3) <- as.matrix(logcounts(sce.sc3))
  counts(sce.sc3) <- as.matrix(counts(sce.sc3))
  
  # NOT repeat !!!! 
  set.seed(2020)
  sce.sc3 <- sc3(sce.sc3, ks = seq(9,19,2), biology = FALSE) # svm_max = 5000 is default!!!
  # Setting SC3 parameters...
  # Your dataset contains more than 2000 cells. Adjusting the nstart parameter of kmeans to 50 for faster performance...
  # Calculating distances between the cells...
  # Performing transformations and calculating eigenvectors...
  # Performing k-means clustering...
  # Calculating consensus matrix...
  
  traceback()
  
  # When the sce.sc3 object is prepared for clustering, SC3 can also estimate the optimal number of clusters k in the dataset
  # NOT repeat, runs 10 mins  !!!!
  sce.sc3 <- sc3_estimate_k(sce.sc3)
  str(metadata(sce.sc3)$sc3)
  # $ k_estimation   : num 19
  
  # to save space, transform back matrix to sparse matrix
  assayNames(sce.sc3)
  #[1]  "counts" "logcounts"
  save(sce.sc3, file='GSE130146_robustness/sce_SC3.RData', compress=TRUE) 
  assays(sce.sc3) <- assays(sce.sc3)[1]
  
  #save(sce.sc3, file='GSE130146_robustness/sce_SC3.RData', compress=TRUE) 
  
  gc()
  # END DO NOT REPET !!!!!!!!!!!!!
  
  ##  writing input fiels for QuanTC -----------------------------------------
  M_9 = (sce.sc3@metadata$sc3$consensus$`9`$consensus)
  M_11 = (sce.sc3@metadata$sc3$consensus$`11`$consensus)
  M_13 = (sce.sc3@metadata$sc3$consensus$`13`$consensus)
  M_15 = (sce.sc3@metadata$sc3$consensus$`15`$consensus)
  M_17 = (sce.sc3@metadata$sc3$consensus$`17`$consensus)
  M_19 = (sce.sc3@metadata$sc3$consensus$`19`$consensus)
  # take average of the Consensus.Cluster-agreeable clustering results to get cell-cell similarity matrix M
  M = (M_9+M_11+M_13+M_15+M_17+M_19)/6 

  QuanTC_input_dir ='F:/projects/QuanTC/QuanTC-modified/Input/mHEP_GSE130146/'
  write.table(M, file=paste0(QuanTC_input_dir,'cell-cell.csv'), 
            row.names=FALSE, col.names = FALSE, sep=',')   

  logmat <- as.matrix(logcounts(sce))
  dim(logmat)
  # [1] 4000 1531
  logmat <- logmat[rownames(sce.sc3),]
  dim(logmat)
  # [1] 3999 1531
  write.table(logmat, file=paste0(QuanTC_input_dir, 'logmat.txt'), 
              sep='\t', row.names=FALSE, col.names=FALSE)
  
  true_label = as.vector(colData(sce)$label)
  write.table(true_label, file=paste0(QuanTC_input_dir, 'true_label.txt'), 
              sep='\t', row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  gene_name = rowData(sce)[rownames(logmat),]$Symbol
  length(gene_name)
  #[1] 3999
  gene_name[1:4]
  #[1] "Mesp1"  "Dppa5a" "Phlda2" "Tdgf1" 
  write.table(gene_name, file=paste0(QuanTC_input_dir, 'gene_name.txt'), 
              sep='\t', row.names=FALSE, col.names = FALSE, quote=FALSE)
  
}


################################################################################  
## section 4) applying cell clustering with different methods                 ##
################################################################################  
#load(file='GSE130146_robustness/sce.GSE130146_noenderdormPgerm.RData')
load(file='../data/sce.GSE130146_noenderdormPgerm.RData')
sce
# class: SingleCellExperiment 
# dim: 4000 1531 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(4000): ENSMUSG00000030544 ENSMUSG00000060461 ...
# ENSMUSG00000070814 ENSMUSG00000031951
# rowData names(2): ID Symbol
# colnames(1531): AAACCTGAGTTTCCTT-1 AAACCTGCATGAAGTA-1 ...
# TTTGTCATCGGCTTGG-1 TTTGTCATCTGCTTGC-1
# colData names(4): Sample Barcode sizeFactor label
# reducedDimNames(3): PCA TSNE UMAP
# altExpNames(0):

# label is the results of SNNGraph_k5 in section 2
sce$C_SNNGraph_k5 = sce$label

  {
    library(dplyr)
    library(scater)
    library(scran)
    library(SC3)
    library(Seurat) #Seurat version 4.0.6
    #library(leiden) 
    #library(monocle3)
    
    
    subDir = '../data/'  #'GSE130146_robustness/'
    QuanTCDir = 'QuanTC_Output/k8/'
    
    {
      parameters = list()
      parameters$k = 10 # An integer number of nearest neighboring cells to use when creating the k nearest neighbor graph for Louvain/Leiden/SNNGraph clustering.
      
      table(colData(sce)$C_SNNGraph_k5)
      #  1  10  11  12  14  17   3   4   5   6   7   8   9 
      # 102  38 208  34  39  11 289  98  86  75 151 160 240 
      
      ########################################################################################
      # 8.1) # extract SNNGraph clusters (by scran) with two settings for the parameter k
      # The parameter k indicates the number of nearest neighbors to consider during graph construction, whicc
      # we set to 5 for small number of cells (e.g., <500) and 10 to large number of sequenced cells (e.g., 4k).
      # https://nbisweden.github.io/single-cell_sib_scilifelab/session-clustering/clustering.html
      ########################################################################################
     { 
        # k: An integer scalar specifying the number of nearest neighbors to consider during graph construction.
        SNNGraph.ID <- scran::buildSNNGraph(sce, k= parameters$k, use.dimred = 'PCA')
        SNNGraph.ID <- igraph::cluster_walktrap(SNNGraph.ID)$membership
        
        # check the agreeement between new and original clusters using the SNNGraph method
        table(as.vector(sce$label), SNNGraph.ID)
        #    SNNGraph.ID
        # 1   2   3   4   5   6   7   8   9
        # 1    0  11   0   0   0   0  89   2   0
        # 10  36   0   0   1   0   0   0   1   0
        # 11   0   0  86 120   0   0   0   2   0
        # 12   0   2   0   0  32   0   0   0   0
        # 14   0   0   0   0   0  39   0   0   0
        # 17   0   0   2   0   0   9   0   0   0
        # 3    1   0   0 284   0   0   0   4   0
        # 4    0   0   4  56   0   5   0  33   0
        # 5    0   0   5   0   0   0   0   3  78
        # 6    3  67   0   0   5   0   0   0   0
        # 7   68  43   0  20   0   0   2  18   0
        # 8    0   0 137   1   0  13   0   0   9
        # 9   13   2   4  34   0   0  13 174   0
        colData(sce)$C_SNNGraph_k10 = factor(SNNGraph.ID)
        
        
        #save(sce, file=paste0(subDir,'/sce.GSE130146_noenderdormPgerm.RData'), compress=TRUE) # !!!!!!!!!!!!!!!!!
      } 
      
      ################################################################
      # 8.2) # extract consensus clusters  
      # refer to http://127.0.0.1:11637/library/SC3/doc/SC3.html
      # needs to customize ks accordingily per dataset, the larger range the longer running time.
      # In this case, ks=3,5,7 are tested, 
      # and for the related soft-thresholding clustering (QuanTC method), 
      # we had take average of the Consensus.Cluster-agreeable clustering results of k=3,5,7 to get cell-cell similarity matrix M
      ##################################################################
      
      load(paste0(subDir,'/sce_SC3.RData') ) 
      {
        sce.sc3  # optimale num =19 by SC3
        
      # load(file='sce_E8.25_HEP.RData')
        colData(sce)$C_consensus_ks13 = colData(sce.sc3)$sc3_13_clusters
        colData(sce)$C_consensus_ks15 = colData(sce.sc3)$sc3_15_clusters
        colData(sce)$C_consensus_ks17 = colData(sce.sc3)$sc3_17_clusters
        colData(sce)$C_consensus_ks19 = colData(sce.sc3)$sc3_19_clusters
        
        rm(sce.sc3)
        
        #save(sce, file=paste0(subDir,'/sce.GSE130146_noenderdormPgerm.RData'), compress=TRUE) # !!!!!!!!!!!!!!!!!
      } 
      
      ################################################################
      # 4.3) # extract Leiden clustering (using Seurat)  
      # https://satijalab.org/seurat/articles/get_started.html
      # Leiden requires the leidenalg python.
      # We apply the Seurat packge, by setting the algorithm =4 in the function FindClusters()	
      # This parameter decides the algorithm for modularity optimization 
      # (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 
      # 3 = SLM algorithm; 4 = Leiden algorithm). 
      # The resolution parameter: use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
      # Seurat author recommended that:
      # We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
      ###################################################################
      # generate a pseudo count to run Seurat
      {
      
        # convert from SingleCellExperiment
        sce.seurat <- as.Seurat(sce)
        # Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
        sce.seurat
        
        # Computes the k.param nearest neighbors for a given dataset
        ndim = dim(sce.seurat[['PCA']])[2]
        sce.seurat <- FindNeighbors(sce.seurat, reduction = "PCA", k.param = parameters$k, dims = 1:ndim)
        
        sce.seurat <- FindClusters(sce.seurat, resolution = 0.4, algorithm = 4) # smaller  number of communities 
        table(Idents(sce.seurat), as.vector(colData(sce)$label))  
        #     1  10  11  12  14  17   3   4   5   6   7   8   9
        # 1   0   2  73   0   0   0 283  51   0   0  21   0  20
        # 2   0   0 129   0   0   0   1   1  73   0   0 131   6
        # 3   6   3   6   0   0   0   3  42  13   0  17   4 201
        # 4   3  32   0   0   0   0   2   0   0  68 109   0  11
        # 5  93   1   0   2   0   0   0   0   0   0   4   0   2
        # 6   0   0   0   0  39  11   0   4   0   0   0  25   0
        # 7   0   0   0  32   0   0   0   0   0   7   0   0   0
        
        colData(sce)$C_Leiden_0.4 = Idents(sce.seurat)
        
        sce.seurat <- FindClusters(sce.seurat, resolution = 0.8, algorithm = 4) # smaller  number of communities 
        table(Idents(sce.seurat), as.vector(colData(sce)$label))  
        #      1  10  11  12  14  17   3   4   5   6   7   8   9
        # 1    0   0 127   0   0   0   2   0   0   0   0  69   6
        # 2    6   0   4   0   0   0   0   7  11   0  11   3 150
        # 3    3   0   0   0   0   0   1   0   0  67  84   0  12
        # 4    0   0  11   0   0   0 132   0   0   0  12   0   4
        # 5    0   0   2   0   0   0   0   1  73   0   0  57   1
        # 6    0   0   9   0   0   0  26  79   2   0   0   1   4
        # 7    0   0  54   0   0   0  44   0   0   0   0   0  16
        # 8   93   1   0   2   0   0   0   0   0   0   6   0   2
        # 9    0   0   0   0  39  11   0   8   0   0   0  29   0
        # 10   0   2   1   0   0   0  75   1   0   0   5   1   1
        # 11   0   5   0   0   0   0   8   2   0   0   6   0  44
        # 12   0  30   0   0   0   0   1   0   0   1  27   0   0
        # 13   0   0   0  32   0   0   0   0   0   7   0   0   0
        
        colData(sce)$C_Leiden_0.8 = Idents(sce.seurat)
        
        sce.seurat <- FindClusters(sce.seurat, resolution = 1.2, algorithm = 4) # smaller  number of communities 
        table(Idents(sce.seurat), as.vector(colData(sce)$label))  
        #      1  10  11  12  14  17   3   4   5   6   7   8   9
        # 1    0   0 128   0   0   0   1   0   0   0   0  61   6
        # 2    6   0   5   0   0   0   0   0  11   0  11   3 142
        # 3    0   0  11   0   0   0 132   0   0   0  12   0   4
        # 4    0   0   0   0   0   0   0   1  72   0   0  65   1
        # 5    0   0  54   0   0   0  44   0   0   0   0   0  16
        # 6    3   0   0   0   0   0   0   0   0  67  37   0   6
        # 7    0  30   0   0   0   0   2   0   0   1  73   0   6
        # 8   93   1   0   2   0   0   0   0   0   0   7   0   2
        # 9    0   0   7   0   0   0  29  50   0   0   0   0   3
        # 10   0   2   1   0   0   0  75   1   0   0   5   1   1
        # 11   0   0   0   0  39   1   0   4   0   0   0  29   0
        # 12   0   5   0   0   0   0   6   1   0   0   6   0  45
        # 13   0   0   2   0   0   0   0  41   3   0   0   1   8
        # 14   0   0   0  32   0   0   0   0   0   7   0   0   0
        # 15   0   0   0   0   0  10   0   0   0   0   0   0   0
        
        colData(sce)$C_Leiden_1.2 = Idents(sce.seurat)
        
       rm(sce.seurat)
        
      } 
      save(sce, file=paste0(subDir,'/sce.GSE130146_noenderdormPgerm.RData'), compress=TRUE)  
      
      ################################################################
      # 4.4 ) # extract the QUANTC-assigned clusters
      # k=11 was optimalized by QuanTC pipeline
      # NOT runnable over night !!!!!!!!!!!
      ##################################################################
      {
        C_TC <- read.table(paste0(QuanTCDir,'C_TC.txt'))
        C_TC <- C_TC[,1]
        length(C_TC) #[1] 1531

        index_TC <- read.table(paste0(QuanTCDir,'index_TC.txt'))
        index_TC <- index_TC[,1]
        unique(C_TC[index_TC]) # 9 verified the C_TC is the cluster ID generated by QuanTC

        ## replace the QuanTC.cluster IDs, TC is the last, to be consistented with those shoing in Fig S2
        tmp <- data.frame(C_TC)
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 1, 'C1'))
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 2, 'C2'))
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 3, 'C3'))
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 4, 'C4'))
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 5, 'C5'))
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 6, 'C6'))
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 7, 'C7'))
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 8, 'C8'))
        tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 9, 'TC'))

        table(as.vector(sce$C_SNNGraph_k5), tmp[,1])
        #     C1  C2  C3  C4  C5  C6  C7  C8  TC
        # 1    0   0  20   0   1   0   0  76   5
        # 10   0   0   0   0   3   0   0   0  35
        # 11  13  54   0  29   4   0   3   0 105
        # 12   0   0  31   0   0   0   0   2   1
        # 14   0   0   0   0   0   0  39   0   0
        # 17   0   0   0   0   0   0   9   0   2
        # 3   96   3   0 126  12   0   1   0  51
        # 4   16   0   0   0  68   0   3   0  11
        # 5    0   4   0   0   4  67   0   1  10
        # 6    0   0  72   0   0   0   0   0   3
        # 7   37   0  29   1   1   0   0   4  79
        # 8    0  88   0   1   8  24  16   0  23
        # 9   68   5   5   2  25   1   1  16 117

        colData(sce)$C_Soft <- tmp[,1]
      }
      
      save(sce, file=paste0(subDir,'/sce.GSE130146_noenderdormPgerm.RData'), compress=TRUE) # !!!!!!!!!!!!!!!!!
      
      # 
      ## 4.5)  plot different clustering restuls                              
      ####################################################################
      library(scater)
      if(!'TSNE' %in% reducedDimNames(sce)) sce <- runTSNE(sce, dimred="PCA")
      
      {  
        x <- grep('C_', colnames(colData(sce)))
        (n=length(x)) # 10
        colnames(colData(sce))[x]
        # [1] "C_SNNGraph_k5"    "C_SNNGraph_k10"   "C_consensus_ks13" "C_consensus_ks15"
        # [5] "C_consensus_ks17" "C_consensus_ks19" "C_Leiden_0.4"     "C_Leiden_0.8"    
        # [9] "C_Leiden_1.2"     "C_Soft"
        # 
        pdf(file="TSNE_clustering_methods.pdf", width=10, height=9)
        gridExtra::grid.arrange(
          plotReducedDim(sce, dimred='TSNE', colour_by='C_SNNGraph_k5', #add_legend=FALSE,
                         text_by='C_SNNGraph_k5', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('C_SNNGraph_k5')  #+ ylim(5,20)  
          ,plotReducedDim(sce, dimred='TSNE',colour_by='C_SNNGraph_k10', #add_legend=FALSE,
                          text_by='C_SNNGraph_k10', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('C_SNNGraph_k10') #+ ylim(5,20)
          ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks13', #add_legend=FALSE,
                          text_by='C_consensus_ks13', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('C_consensus_ks13')  #+ ylim(5,20)
          ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks15', #add_legend=FALSE,
                          text_by='C_consensus_ks15', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('C_consensus_ks15')  #+ ylim(5,20)
          ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks17', #add_legend=FALSE,
                          text_by='C_consensus_ks17', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('consensus_ks17') #+ ylim(5,20) 
          ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks19', #add_legend=FALSE,
                          text_by='C_consensus_ks19', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('consensus_ks19') #+ ylim(5,20)
          ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Leiden_0.4', #add_legend=FALSE,
                          text_by='C_Leiden_0.4', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('Leiden_0.4') #+ ylim(5,20)
          ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Leiden_0.8', #add_legend=FALSE,
                          text_by='C_Leiden_0.8', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('Leiden_0.8') #+ ylim(5,20)
          ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Leiden_1.2', #add_legend=FALSE,
                          text_by='C_Leiden_1.2', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('Leiden_1.2') #+ ylim(5,20)
          ,ncol=3
        )
        gridExtra::grid.arrange(
          plotReducedDim(sce, dimred='TSNE', colour_by='C_Soft', #add_legend=FALSE,
                         text_by='C_Soft', text_size = 4, text_colour='black', point_size=0.5) + 
            ggtitle('C_Soft')  #+ ylim(5,20)  
          ,ncol=3, nrow=3)
        
        dev.off()
        
      }  
    }
    
    
}
