## This code calculates consensus matrixs to run QuanTC
## last upded 8/10/2021
## by Holly Yang

# R4.0.2
setwd('F:/projects/BioTIP/result/QuanTC_simulation/')


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
 library(bluster)

 library(ggplot2)
 library(gridExtra)  # required by grid.arrange()
 library(ggpubr)  # required by p-values on bar
 
 library(SC3)
 packageVersion('SC3') # 1.18.0
 
 library(MouseGastrulationData)
 head(EmbryoCelltypeColours)
 length(EmbryoCelltypeColours) # 37
 
 inputpath = "Input/" 
 BioTIPoutputpath = "BioTIP/"
 ICoutputpath = "original_IC/"
 
##############################################################
## 1) translate the MATHLAB data into an R object           ##
## original code: 2_QuanTC_simulation_cluster_notused.R     ##
            ## DO NOT repeat    
##############################################################
# refer to 1_export_data.m
 X <- NULL # a matrix of 18 gene x 5363 cells
 for(i in 1:5){
   fileName = paste0(inputpath,'data',i,'.txt')
   tmp <- read.table(fileName, sep=",")
   dim(tmp) # 312 18,  a cell x gene matrix
   X <- cbind(X, t(tmp))
 } 
 
 # check normalitity
 # boxplot(X)
 # boxplot(log2(X+1))
 hist(log2(X+1), 50, ylab='log2(x+1)', main="data of 5363 cells")
 dev.copy2pdf(file='hist_data_18gene_exp.pdf')
 
 y <- NULL # a matrix of 18 gene x 5363 cells
 for(i in 1:5){
   fileName = paste0(inputpath,'true_label',i,'.txt')
   tmp <- read.table(fileName, sep=",")
   dim(tmp) # 312 1 
   y <- c(y, tmp[,1])
 } 
 length(y)
 # [1] 5363
 table(y)
 #   1    2    3    4    5 
 # 577 1849  801  985 1151 
 tmp = y
 tmp[tmp==1] = 'E'
 tmp[tmp==2] = 'I1'
 tmp[tmp==3] = 'I2'
 tmp[tmp==4] = 'M'
 tmp[tmp==5] = 'TC'
 tmp = factor(tmp, levels = c('E','I1', 'I2', 'M', 'TC'))
 y = factor(y, levels=c('1','2','3','4','5'))
 
 fileName = paste0(inputpath,'gene_name.txt')
 gene_name <- read.table(fileName, sep=",")
 gene_name=as.character(gene_name[1,])
 
 rownames(X) <- gene_name
 colnames(X) <- y
 
 sce <- SingleCellExperiment(list(counts=X, logcounts=log2(X+1)),
                             colData=DataFrame(true_label=y,
                                               cell_type=tmp),
                             rowData=DataFrame(gene_name=gene_name),
                             metadata=list(study="EMT simulation genearted by QuanTC authors")
 )
 
save(sce, file='sce.RData')
 
######################################################################## 
 load('sce.RData')
 sce
 # class: SingleCellExperiment 
 # dim: 18 5363 
 # metadata(1): study
 # assays(2): counts logcounts
 # rownames(18): snailt SNAIL ... Ncad Ovol2
 # rowData names(1): gene_name
 # colnames(5363): 1 2 ... 5362 5363
 # colData names(3): true_label cell_type C_SNNGraph.k200
 # reducedDimNames(4): PCA PCA.rawcount TSNE UMAP
 # altExpNames(0):
     
  
 
 
##################################################
## 2) Estimate cell-cell similarity by SC3      ##
# http://127.0.0.1:11637/library/SC3/doc/SC3.html
##################################################
 
 table(colLabels(sce))  
 #   1    2    3    4 
 # 2373 1017  752 1221 
 
 
 # define feature names in feature_symbol column
 rowData(sce)$feature_symbol <- rownames(sce)
 # remove features with duplicated names
 any(duplicated(rowData(sce)$feature_symbol))  # F
 range(counts(sce))
# [1]      3.927686e-20 6.854571e+02
 
 ### to run SC3 successfully, transform sparsematrix to matrix  !!!!!!!!!!!!!
 counts(sce) <- as.matrix(counts(sce))
 logcounts(sce) <- as.matrix(logcounts(sce))
 
 sum(logcounts(sce)<1e-16,2)/nrow(logcounts(sce))>min_exp
 
 
 
# NOT repeat, runs 1 hour !!!! 
# biology: boolean parameter, defines whether to compute differentially expressed genes, marker genes and cell outliers.
 set.seed(2020)
 sce <- sc3(sce, ks = 3:6, biology = TRUE, gene_filter = FALSE, svm_max = 6000)  
 # Setting SC3 parameters...
 # Your dataset contains more than 2000 cells. Adjusting the nstart parameter of kmeans to 50 for faster performance...
 # Calculating distances between the cells...
 # Performing transformations and calculating eigenvectors...
 # Performing k-means clustering...
 # Calculating consensus matrix...
 # Calculating biology...
 
 traceback()
 
 # When the sce object is prepared for clustering, SC3 can also estimate the optimal number of clusters k in the dataset
 # NOT repeat, runs 10 mins  !!!!  
 sce <- sc3_estimate_k(sce)
 str(metadata(sce)$sc3)
 # $ k_estimation   : num 1
 
 sce
 #: SingleCellExperiment 
 # dim: 18 5363 
 # metadata(2): study sc3
 # assays(2): counts logcounts
 # rownames(18): snailt SNAIL ... Ncad Ovol2
 # rowData names(19): gene_name feature_symbol ... sc3_5_de_padj
 # sc3_6_de_padj
 # colnames(5363): 1 2 ... 5362 5363
 # colData names(11): true_label cell_type ... sc3_5_log2_outlier_score
 # sc3_6_log2_outlier_score
 # reducedDimNames(4): PCA PCA.rawcount TSNE UMAP
 # altExpNames(0):    

 # to save space, transform back matrix to sparse matrix
 
 counts(sce) <- as(counts(sce), 'dgCMatrix')
 logcounts(sce) <- as(logcounts(sce), 'dgCMatrix')

 save(sce, file='sce_SC3.RData')    ##  save again !!!!!!!!!!!!!!!!!!!
 
 col_data <- colData(sce)
 head(col_data[ , grep("sc3_", colnames(col_data))], 3)
 
 row_data <- rowData(sce)
 head(row_data[ , grep("sc3_", colnames(row_data))], 3)
 
 fileNames <- c('PCA','TSNE', 'UMAP')
 reducedDimNames <- c('PCA', 'TSNE', 'UMAP')

 for(i in 1:length(fileNames)){
    plist=list()
    plist[[1]] <-  plotReducedDim(sce, reducedDimNames[i], colour_by="label",
                                  text_by='label', text_size = 5,
                                  point_size=1)
    plist[[2]] <- plotReducedDim(sce, reducedDimNames[i], colour_by = "sc3_3_clusters", 
                                 text_by='sc3_3_clusters', text_size = 5,
                                 point_size=1 ) #+ #"size_by =sc3_3_log2_outlier_score") +
      # scale_color_manual(values= as.vector(EmbryoCelltypeColours))
    plist[[3]] <- plotReducedDim(sce, reducedDimNames[i], colour_by = "sc3_4_clusters", 
                                 text_by='sc3_4_clusters', text_size = 5, 
                                 point_size=1 ) #+ #size_by ="sc3_3_log2_outlier_score") + 
      # scale_color_manual(values= as.vector(EmbryoCelltypeColours))
    plist[[4]] <- plotReducedDim(sce, reducedDimNames[i], colour_by = "sc3_5_clusters", 
                                 text_by='sc3_5_clusters', text_size = 5, 
                                 point_size=1 ) #+ #size_by ="sc3_3_log2_outlier_score") + 
      # scale_color_manual(values= as.vector(EmbryoCelltypeColours))
    plist[[5]] <- plotReducedDim(sce, reducedDimNames[i], colour_by = "sc3_6_clusters",
                                 text_by='sc3_6_clusters', text_size = 5, 
                                 point_size=1 ) #+ #size_by ="sc3_3_log2_outlier_score") + 
      # scale_color_manual(values= as.vector(EmbryoCelltypeColours))
    
    pdf(file=paste0('SC3_kselection_',fileNames[i],'.pdf'), width=14, height=9)
    do.call("grid.arrange", c(plist, ncol=3))
    dev.off()
    
 } 
 
 
 M_3 = (sce@metadata$sc3$consensus$`3`$consensus)
 M_4 = (sce@metadata$sc3$consensus$`4`$consensus)
 M_5 = (sce@metadata$sc3$consensus$`5`$consensus)
 M_6 = (sce@metadata$sc3$consensus$`6`$consensus)
 
 
 # take average of the Consensus.Cluster-agreeable clustering results to get cell-cell similarity matrix M
 # because that k=20 is the best matching by eye-ball
 gc()
 M = (M_3+M_4+M_5+M_6)/4   
 write.table(M, file='QuanTC_Input_M_cell-cell.csv', row.names=FALSE, col.names=FALSE, sep=',') # !!!!!!!!!!!!!!!!!!!!!!!!!
 system("head QuanTC_Input_M_cell-cell.csv")

 
pdf(file='sc3_plot_markers.k4.pdf', height=20)
sc3_plot_markers(sce, k = 4,
                   show_pdata = c(
                     "label" ,
                     "sc3_4_clusters"
                   )
 )
 dev.off()
 table(col_data$cell_type, col_data$sc3_4_clusters)
 #       1    2    3    4
 # E    35    6    0  536
 # I1 1848    0    0    1
 # I2  110    0    0  691
 # M    55    0  927    3
 # TC 1088    1    2   60
 
 
 
 pdf(file='sc3_plot_markers.k5.pdf', height=20)
 sc3_plot_markers(sce, k = 5,
                  show_pdata = c(
                      "label" ,
                      "sc3_5_clusters"
                  )
 )
 dev.off()
 table(col_data$cell_type, col_data$sc3_5_clusters)
 #    1    2    3    4    5
 # E     7   33    0    0  537
 # I1    0 1843    0    6    0
 # I2    0   59    1  741    0
 # M     0   50  935    0    0
 # TC    1 1026   60   18   46
 
 
 
 pdf(file='sc3_plot_markers.k6.pdf', height=20)
 sc3_plot_markers(sce, k = 6,
                  show_pdata = c(
                      "label" ,
                      "sc3_6_clusters"
                  )
 )
 dev.off()
 table(col_data$label, col_data$sc3_6_clusters)
 #       1    2    3    4    5    6
 # E     8    0    4    0  533   32
 # I1    0 1768   46    0    0   35
 # I2    0    0   38    0    0  763
 # M     0    0    7  927    0   51
 # TC    1    5  560    2   36  547
 
 
 
 
 