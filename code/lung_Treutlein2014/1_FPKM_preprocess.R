#### README ############
## There are 9 sections in this code 
## Section 1) download and process the FPKM and sample infor from GEO            
## section 2) Quanlity control and feature selection                  
## section 3) clustering cells by gene markers          
## section 4) clustering cells without gene markers based on 10.3k expressed genes that were used to construct trajectory by Monocle3        
## Section 5 (OPTINAL): constructing trajectory using Monocle 
## Section 6) select Global ~3k HVG  
## Section 7) Prepare inputs for QuanTC with 3198 genes and 131 cells           
## section 8) cluster 131 cells of 4531 genes using different methods, and plot the cell clusters              
## scetion 9) construct and visualize teh trajectory using scater package
## last update 2/8/2022
## by Holly Yang  ##############

setwd('F:/projects/BioTIP/result/GSE52583')

####################################################################################
### Section 1 ) download and process the FPKM and sample infor from GEO           ##
## columns "age" cells" "cellName" are from GEO                                   ##
## column "putative_cell_type" are from the publishgd Stable 3                    ##
## factor "cell_type" merges both age and putative_cell_type                      ##
## "CellType" is added in Section 4, the biomarker-based cell clusters            ##
## "Cluster" is added in Section 5, the unsupervised cell clusters (k=6)          ##
####################################################################################
{
  library(GEOquery)
  GSE52583 <- getGEO(GEO = 'GSE52583', filename = NULL, destdir = '../../data/GSE52583_lung_cellfate_usedbyMojtahedi2016/GSE52583GEOquery',
                     GSElimits = NULL, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = TRUE)
  class(GSE52583)  #[1] "list"
  names(GSE52583)  
  #  [1] "GSE52583-GPL13112_series_matrix.txt.gz" "GSE52583-GPL16417_series_matrix.txt.gz"
  class(GSE52583[[1]])  #[1] "ExpressionSet"
  save(GSE52583, file= '../../data/GSE52583_lung_cellfate_usedbyMojtahedi2016/GSE52583GEOquery/FPKM.GSE52583_ExpressionSet_list.rData',compress=T)
  
  ## generate meta_table for cells ################
  cli <- rbind(pData(GSE52583[[1]]), pData(GSE52583[[2]]))
  dim(cli)     # [1] 201  46
  toDelete <- NULL
  for(i in 1:ncol(cli)) {
    if(length(table(cli[,i])) ==1) toDelete <- c(toDelete, i) 
  }
  cli <- cli[,-toDelete]
  dim(cli)  #[1] 201  10
  cli$age <- unlist(lapply(cli$characteristics_ch1.2, function(x) unlist(strsplit(as.character(x), split=" day "))[2]))
  cli$age[which(cli$age=='107')]='Adult'
  cli$genotype <- unlist(lapply(cli$characteristics_ch1.3, function(x) unlist(strsplit(as.character(x), split=": "))[2]))
  table(cli$age, cli$genotype )
  #        Sftpc-Cre-ERT2-rtta -/- tetO-HIST1H2BJ-GFP+/-) wild type
  # 14.5                                               0        45
  # 16.5                                               0        27
  # 18.5                                               0        83
  # Adult                                             46         0
  
  cli$replicate <- unlist(lapply(cli$title, function(x) unlist(strsplit(as.character(x), split=", "))[2]))
  cli$cells <- unlist(lapply(cli$title, function(x) unlist(strsplit(as.character(x), split=", "))[3]))
  table(cli$cells)
  #200 cell bulk control       no cell control          single cells 
  #                    2                     1                   198 
  rownames(cli) <- cli$geo_accession
  cli <- cli[,c]
  tmp <- unlist(lapply(cli$supplementary_file_1, function(x) unlist(strsplit(as.character(x), split="suppl/"))[2]))
  tmp <- unlist(lapply(tmp, function(x) unlist(strsplit(x, split="_IL"))[1]))
  tmp1 <- unlist(lapply(tmp, function(x) unlist(strsplit(x, split="_"))[2]))
  tmp2 <- unlist(lapply(tmp, function(x) unlist(strsplit(x, split="_"))[3]))
  tmp3 <- unlist(lapply(tmp, function(x) unlist(strsplit(x, split="_"))[4]))
  
  cli$cellName <- paste(tmp1,tmp2,tmp3,sep="_")
  colnames(cli)[8] <- 'biosample'
  colnames(cli)[9] <- 'SRX'
  
  cli <- cli[,c(6,7,11:15)]
  
  ##  Treutlein2014 suppl data 3,  match to cli  
  Data3 <- read.table('../../data/GSE52583_lung_cellfate_usedbyMojtahedi2016/Treutlein2014/Treurlein2014_SData3.txt',header=T,sep='\t')
  dim(Data3)  #[1]    82 23275
  colnames(Data3)[1:10]
  # [1] "cell_name"          "time_point"         "sample"            
  # [4] "putative_cell_type" "X0610005C13Rik"     "X0610007C21Rik"    
  # [7] "X0610007L01Rik"     "X0610007N19Rik"     "X0610007P08Rik"    
  # [10] "X0610007P14Rik"
  table(Data3$putative_cell_type)
  #     AT1      AT2       BP     bulk    ciliated    Clara 
  #      41       12       13        2        3       11 
  
  cli$putative_cell_type <- Data3[match(cli$cellName,Data3$cell_name),]$putative_cell_type
#  write.table(cli, file='GSE52583_sample_annotation_xy.txt', sep='\t')
  
  ## generate matrix of FPKM ######################
  myDir <- "F:/projects/BioTIP/data/GSE52583_lung_cellfate_usedbyMojtahedi2016/GSE52583_RAW/" 
  COLs <- c('tracking_id', 'class_code', 'nearest_ref_id', 'gene_id', 'gene_short_name', 'tss_id', 
          'locus', 'length','coverage','FPKM','FPKM_conf_lo','FPKM_conf_hi','FPKM_status')

  files <- list.files(path = myDir, pattern='*.gz')
  (n <- length(files))  # 201
  tmp <- read.table(file=paste0(myDir,files[1]),sep='\t',header=FALSE, comment="")
  colnames(tmp) <- COLs

  FPKM <- NULL
  for(i in 1:n)
  {
    tmp <- read.table(file=paste0(myDir,files[i]),sep='\t',header=FALSE, comment="")
    colnames(tmp) <- COLs
    FPKM <- cbind(FPKM, tmp$FPKM)
  }
  rownames(FPKM) <- tmp$gene_id
  colnames(FPKM) <-  unlist(lapply(files, function(x) unlist(strsplit(x, split="_"))[1]))
  dim(FPKM)   #[1] 23837   201
  
  write(rownames(FPKM), file= '../../data/GSE52583_lung_cellfate_usedbyMojtahedi2016/Treutlein2014/FPKM_correctSymbol.txt')
#  save(FPKM, file='FPKM_matrix_nofilter.rData', compress=T)
  any( rownames(cli)!=colnames(FPKM))  #[1] FALSE
  
}

################################################################################################################
## section 2) Quanlity control and feature selection                  
## GEO downloaded FPKM matrix                                                                       23837   201
## remove_spike, mito, non-annotated transcripts                                                    23231   201
## collapse the FPKM values for duplicated symbol by the mean value                                 22854   201
## cell_quality_filtering by removing 
#### either very low mRNA recovery or potential doublets or triplets                                22854   196
## feature selection based on mean-variance relationship in trajectory construction                 10359   196
## feature selection (in section 7) of coding gene, miRNA and lincRNA  (expressed)                  10251   196
#### feature selection (in section 7) of HVGs                                                       3198   196
## focusing on cells along the AT2 trajectory                                                       3198   131
##################################################################################################################
library(monocle)
{
  cell_median_FPKM <- apply(FPKM, 2, function(df) median(df, na.rm=TRUE) )
  table(cell_median_FPKM==0) # TRUE 201 !!!!  FPKM has been centralized per cell !!!
  
  genes.ERCC <- grep("ERCC", rownames(FPKM),value=TRUE)
  length(genes.ERCC) # 92
  FPKM.ERCC <- FPKM[genes.ERCC,]
  dim(FPKM.ERCC)   # 92  201
  save(FPKM.ERCC, file='FPKM.ERCC.RData') #!!!!!!!!!!!!
  
  ### housekeeping genes given in extend Data Fig 3
  HK <- c("Gusb", "Tbp","Ppih", "Tfrc", "Sdha", "Pgk1", "B2m", "Ldha", "Gapdh", "Hsp90ab1", "Actb") 
  all(HK %in% rownames(FPKM))  #[1] TRUE
  
  
## 2.1) prepare annotate table ------------------------------
{  # https://ivanek.github.io/analysisOfGenomicsDataWithR/03_AnnotationResources_html.html
  # https://www.biostars.org/p/147351/
  # https://www.biostars.org/p/147351/
  library(biomaRt)   #biomaRt_2.42.0
  ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
  grep("Synonyms",listAttributes(ensembl), value=T)
  annot<-getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), 
               filters= "mgi_symbol", value=rownames(FPKM), mart=ensembl)
  dim(annot)  #[1] 20990     7
  table(rownames(FPKM) %in% annot$mgi_symbol)
  # FALSE  TRUE 
  # 2875 20962
  x <- which(!rownames(FPKM) %in% annot$mgi_symbol); length(x) # 2875
  # load the annotation database
  # set up your query genes
  queryGeneNames <- rownames(FPKM)[x]
  
  # use sql to get alias table and gene_info table (contains the symbols)
  # first open the database connection
  #library(DBI)    #DBI_1.1.0
  dbCon <- org.Mm.eg_dbconn()
  # write your SQL query
  sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
  # execute the query on the database
  aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
  dim(aliasSymbol)  #[1] 149181      5
  
  # subset to get your results
  result <- aliasSymbol[which(aliasSymbol[,2] %in% queryGeneNames),c(2,5)]
  dim(result)   # [1] 2476    2
  length(unique(result[,1]) )  #[1] 2430
  unique(result[duplicated(result[,1]),1])
  #[1] "Ang3"          "Egfbp2"        "Aim1"          "Stra13"        "Odz2"         
  # [6] "Odz1"          "Odz3"          "Prl2c4"        "Nat6"          "Rab1"         
  #[11] "Sip1"          "1700058G18Rik" "Abp1"          "Klra22"        "B3gnt1"       
  #[16] "G6b"           "ORF61"         "Dbc1"          "Clca4"         "Siglec5"      
  #[21] "6430411K18Rik" "Mll2"          "6430706D22Rik" "A730008H23Rik" "2810408M09Rik"
  #[26] "Stxbp3a"       "Nrp"           "Dear1"         "Duxbl"         "Plac9"        
  #[31] "H2afb1"        "Gmcl1l"        "Cldn25"        "Cml3"          "Mir193"       
  #[36] "5430421N21Rik" "Snord116"      "Rbmy1a1"    
  # manually correct unrecognized symbles with alias Symbol #!!!!!!!!!!!!!!!!!!!!!!!
  tmp <- unique(c(549,11855,21429,24529,12342,21477, 25735,25741,25746,26420,
                  10086, 35195,38777,35170,4730,11893,68927, 66645, 58335, 80156, 82624,
                  83976, 96333, 56001,59911, 77808,57346, 98574, 107119,42860,107291,
                  116138,116186,116234,116825,125298,125445, 118088, 118604,124424, 
                  119423,12090, 38360, 125467, 125484, 125487))
  result <- result[-which(rownames(result) %in% as.character(tmp)),]
  dim(result) # [1] 2430    2
  length(unique(result[,1]) )  #[1] 2430
  annot2 <- getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), 
                  filters= "mgi_symbol", value=result[,"symbol"], mart=ensembl)  
  annot <- rbind(annot, annot2)
  dim(annot)    # [1] 23213     7
  rm(annot2)
  
#  save(annot, file="ensembl_annot_GSE52583_full.RData", compress=T)
  
  GSE52583.gene.id <- annot$mgi_symbol
  tmp <- match(result$symbol, GSE52583.gene.id)
  length(tmp)   # 2430
  GSE52583.gene.id[tmp[!is.na(tmp)]] <- result$alias_symbol[!is.na(tmp)]
  annot$GSE52583.gene.id <- GSE52583.gene.id
  
  # correct the wrong chromosome_name 
  x <- grep("CHR", annot$chromosome_name); length(x) # 311
  for(i in 1:length(x))
  {
    annot[x[i],"chromosome_name"] <- getBM("chromosome_name", 
                                           filters= "mgi_symbol", value=annot$mgi_symbol[x[i]], mart=ensembl)
  }  
  dim(annot)  # 23213     8
  table(annot$chromosome_name)
  x <- grep("JH58", annot$chromosome_name); length(x) # 6
  for(i in 1:length(x))
  {
    annot[x[i],"chromosome_name"] <- getBM("chromosome_name", 
                                           filters= "mgi_symbol", value=annot$mgi_symbol[x[i]], mart=ensembl)
  }  
  dim(annot)  # 23213     8
  table(annot$chromosome_name)
  x <- grep("GL45", annot$chromosome_name); length(x) # 8
  for(i in 1:length(x))
  {
    annot[x[i],"chromosome_name"] <- getBM("chromosome_name", 
                                           filters= "mgi_symbol", value=annot$mgi_symbol[x[i]], mart=ensembl)
  }  
  dim(annot)  # 23213     8
  table(annot$chromosome_name)
  x <- grep("GL45", annot$chromosome_name); length(x) # 4
  annot[x,"chromosome_name"] <- c(1, "X", "X","X")
  x <- grep("JH58", annot$chromosome_name); length(x) # 3
  annot[x,"chromosome_name"] <- c(4, 5, 4)
  annot[x[1] ,"gene_biotype"] <- "lncRNA"
  
  
  colnames(annot)[which(colnames(annot)=="mgi_symbol")] <- 'gene_short_name'
  annot$strand[which(annot$strand=="1")] ="+"
  annot$strand[which(annot$strand=="-1")] ="-"
  annot$locus <- paste0("chr",annot$chromosome_name,":",annot$start_position,"-",annot$end_position,":",annot$strand) 
  GRanges(annot$locus)
  #GRanges object with 23213 ranges and 0 metadata columns:
#  save(annot, file="ensembl_annot_GSE52583_full.RData", compress=T) # !!!!!!!!!
}  

## 2.2) Remove non-annotated transcripts, including spike, mito transcripts   ------------------------------
{
  FPKM <- FPKM[which(rownames(FPKM) %in% GSE52583.gene.id),]
  dim(FPKM)  #23231    /   filtered: [1] 15897   201
  
 
  ### collapse the FPKM values for duplicated gene symbol by the mean value --------------------------- 
  x <- which(duplicated(rownames(FPKM)))  ; length(x)  #[1] 377
  tmp <- FPKM[x,]
  FPKM <- FPKM[-x,]
  for(i in 1:length(x))
  {
    y <- which(rownames(FPKM)==x[i])
    FPKM[y,] <- apply(rbind(FPKM[y,], tmp[i,]), 2, mean)
  }
  dim(FPKM)  #  22854 /   filtered: [1] 15782   201
}  

## 2.3) arbitarily take the first annotation for the duplicated annotations --------------------
{
  annot <- annot[-which(duplicated(annot$GSE52583.gene.id)),]
  rownames(annot) <- annot$GSE52583.gene.id
  annot <- annot[rownames(FPKM),]
  
  
  pd <- new("AnnotatedDataFrame", data = cli)
  fd <- new("AnnotatedDataFrame", data = annot)
  annot.FPKM <- newCellDataSet(FPKM, phenoData = pd, featureData = fd,
                               lowerDetectionLimit = 1,
                               expressionFamily=tobit())   # Tobits are truncated normal distributions. 

  cds <- relative2abs(annot.FPKM, method = "num_genes") 
  annot.cds <- newCellDataSet(cds, phenoData = pd, featureData = fd,
                              lowerDetectionLimit = 1,
                              expressionFamily=negbinomial.size())   # Negative binomial distribution with fixed variance (which is automatically calculated by Monocle). Recommended
 
 # save(annot.cds, file="GSE52583_annot.cds_allgene.RData", compress=TRUE)  #++++++++++++++++++++++++++++++
}  
  
## 2.4) cell quality control --------------------
  library(monocle)
{  
  valid_cells <- row.names(subset(pData(annot.cds),
                                  cells == "single cells" 
  ))
  annot.cds <- annot.cds[,valid_cells]
  annot.cds <- estimateSizeFactors(annot.cds)
  annot.cds <- estimateDispersions(annot.cds)
  #Removing 291 outliers 
  dim(annot.cds)
  #Features  Samples 
  #   22854      198

  # cut left tail by setting an expression threshold of  2^(-18) !!!!!
  annot.cds <- detectGenes(annot.cds, min_expr = 2^(-18))
  print(head(fData(annot.cds)))
  
  expressed_genes <- row.names(subset(fData(annot.cds),
                                      num_cells_expressed >= 10))
  length(expressed_genes)    #[1] 11333  genes expressed in at least 10 cells

  ## It's also good to look at the distribution of mRNA totals across the cells:
  pData(annot.cds)$Total_mRNAs <- Matrix::colSums(exprs(annot.cds))  #!!!
  summary(pData(annot.cds)$Total_mRNAs)
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  #    728    7703   11310   13023   16333   30088
  
  # We've gone ahead and removed the few cells with either very low mRNA recovery or far more mRNA that the typical cell. ------------------------- 
  # Often, doublets or triplets have roughly twice the mRNA recovered as true single cells, 
  annot.cds <- annot.cds[,pData(annot.cds)$Total_mRNAs < 30000] # !!!!!!!!!!!!!!!
  table(pData(annot.cds)$age)   ## only remove 1 E14.5 cell !!!!!!!!!!!!!!
  # 14.5  16.5  18.5 Adult 
  #   44    27    80    46
  upper_bound <- lower_bound <- NULL
  x <- levels(pData(annot.cds)$age)
  for(i in x)
  {
    tmp <- which(pData(annot.cds)$age==i)
    u_bound <- 10^(mean(log10(pData(annot.cds)$Total_mRNAs[tmp])) +
                     2*sd(log10(pData(annot.cds)$Total_mRNAs[tmp])))
    l_bound <- 10^(mean(log10(pData(annot.cds)$Total_mRNAs[tmp])) -
                     2*sd(log10(pData(annot.cds)$Total_mRNAs[tmp]))) 
    upper_bound <- max(upper_bound,u_bound)           
    lower_bound <- min(lower_bound,l_bound)
  }
  lower_bound  #[1] 1633.736
  upper_bound # [1] 32035.42
  
  qplot(Total_mRNAs, data = pData(annot.cds), color = age, geom =
          "density", ylab="Density") +
    geom_vline(xintercept = lower_bound) +
    geom_vline(xintercept = upper_bound)
  dev.copy2pdf(file="density_Total_mRNAs.pdf")
  
  # so the latter filter is another means of excluding all but single cells from the analysis. 
  # Such filtering is handy if your protocol doesn't allow directly visualization of cell after they've been captured. 
  # Note that these thresholds of 3095 and 38890 mRNAs are specific to this dataset. 
  annot.cds <- annot.cds[,pData(annot.cds)$Total_mRNAs > lower_bound &
                           pData(annot.cds)$Total_mRNAs < upper_bound]
  annot.cds <- detectGenes(annot.cds, min_expr = 0.1)
  dim(annot.cds)
  # Features  Samples 
  #   22854      196 
  table(pData(annot.cds)$age)   ## here remove 2 Adult cells !!!!!!!!!!!!!!
  # 14.5  16.5  18.5 Adult 
  #  45    27    80    44 
  # Log-transform each value in the expression matrix.
  L <- log(exprs(annot.cds[expressed_genes,]))
  any(L==-Inf)  #[1] TRUE
  
  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily
  melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
  
  # verify the data follows a distribution that is roughly lognormal
  # Plot the distribution of the standardized gene expression values.
  qplot(value, geom = "density", data = melted_dens_df) +
    stat_function(fun = dnorm, size = 0.5, color = 'red') +
    xlab("Standardized log(FPKM), 196 single cells") +
    ylab("Density")
  dev.copy2pdf(file="density_Standardized_logcds.pdf")
  
  fData(annot.cds)$expressed_genes <- (row.names(fData(annot.cds)) %in% expressed_genes)
#  save(annot.cds, file="GSE52583_annot.cds_Standardized.RData", compress=TRUE) #++++++++++++++++++++++++++++++
  
  }
  
  ## 2.5) feature selection based on mean gene expression abundance  ----------------------
  ## focus on lncRNA, coding genes, and miRNAs
{
  length(fData(annot.cds)$expressed_genes)  #  11333
  # expressed_genes was only used to exclude two Adult cells  !!!!!!!!!!!!!!!!!!
  # annot.cds <- annot.cds[fData(annot.cds)$expressed_genes, ]
  # dim(annot.cds)
  #Features  Samples 
  #   11333      196 
  table(pData(annot.cds)$putative_cell_type)    
  #     AT1      AT2       BP     bulk ciliated    Clara 
  #      41       12       13        0        3       11  
  pData(annot.cds)$cell_type <- paste(pData(annot.cds)$age, pData(annot.cds)$putative_cell_type, sep="_")  #!!!!
  pData(annot.cds)$cell_type[which(pData(annot.cds)$cell_type=="Adult_NA")] = "Adult_AT2"   #!!!!
  
  pData(annot.cds)$cell_type <- factor(pData(annot.cds)$cell_type, 
                                       levels =c("14.5_NA", "16.5_NA","18.5_BP","18.5_ciliated","18.5_Clara","18.5_AT1","18.5_AT2","Adult_AT2"))
  
  ## Clustering cells (without marker genes) --------------------------
  disp_table <- dispersionTable(annot.cds) # Retrieve a table of values specifying the mean-variance relationship
  dim(disp_table)  #[1] 10365     4
  hist(log(disp_table$mean_expression), 100)
  
  # filter genes based on average expression level --------------------------------
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.01)
  dim(unsup_clustering_genes)    # [1]  10359     4   ######
  
  annot.cds <- setOrderingFilter(annot.cds, unsup_clustering_genes$gene_id)
  table(fData(annot.cds)$use_for_ordering)
  #FALSE  TRUE 
  #12676 10359
  plot_ordering_genes(annot.cds) 
  dev.copy2pdf(file="scarePlot_averageHighgene_Standardized.pdf")
  
 }  
}

###################################################
## section 3) clustering cells by gene markers   ##       
###################################################
{
  AT1_id <- row.names(subset(fData(annot.cds), gene_short_name == "Ager")) # %in% c("Pdpn","Ager","Aqp5") ))
  AT1_id2 <- row.names(subset(fData(annot.cds), gene_short_name == "S100a6")) # %in% c("Pdpn","Ager","Aqp5") ))
  AT1_id3 <- row.names(subset(fData(annot.cds), gene_short_name == "Pdpn")) # %in% c("Pdpn","Ager","Aqp5") ))
  
  AT2_id <- row.names(subset(fData(annot.cds), gene_short_name == "Sftpc"))  #%in% c("Sftpc","Muc1","Sftpb","Abca3","Lyz2")))
  AT2_id2 <- row.names(subset(fData(annot.cds), gene_short_name == "Lyz2"))  #%in% c("Sftpc","Muc1","Sftpb","Abca3","Lyz2")))
  
  cth <- newCellTypeHierarchy()
  cth <- addCellType(cth, "AT1", classify_func =
                       function(x) { x[AT1_id,] > 1 & x[AT1_id2,] > 1 &  x[AT2_id2,] < 1})
  cth <- addCellType(cth, "AT2", classify_func = function(x)
  {  x[AT1_id3,] < 1 & x[AT2_id,] > 1  })
  
  annot.cds <- classifyCells(annot.cds, cth, frequency_thresh =NULL) 
  colnames(pData(annot.cds))
  # [1] "platform_id"         "instrument_model"    "age"                 "genotype"           
  # [5] "replicate"           "cells"               "cellName"            "putative_cell_type" 
  # [9] "remove_by_RECC"      "Size_Factor"         "num_genes_expressed" "Total_mRNAs"        
  # [13] "cell_type"           "CellType"  
  
  table(pData(annot.cds)$CellType, pData(annot.cds)$putative_cell_type)  #!!!!!!!!!!!!
  #            AT1 AT2 BP bulk ciliated Clara
  #  Ambiguous   0   0  1    0        2     1
  #  AT1        38   0  9    0        0     3
  #  AT2         0  11  1    0        1     2
  #  Unknown     3   1  2    0        0     5
}

########################################################
## section 4) clustering cells without gene markers   ##       
########################################################
{
 
  pc <-  plot_pc_variance_explained(annot.cds, return_all = T) # norm_method='log'
  plot(pc$p)
  pc$variance_explained[1:3]*100
  # [1]   8.703065 3.243314 2.495040
  
  annot.cds <- reduceDimension(annot.cds, max_components = 3, num_dim = 6,
                               norm_method = "log",  # because the exprs values are thedownloaded FPKM, log transform is applied  !!!!!!!!!!!!!!!!!!!
                               reduction_method = 'tSNE', verbose = T)
  colnames(pData(annot.cds))  # NO change
  
  annot.cds <- clusterCells(annot.cds, num_clusters = 4)   

#   save(annot.cds, file="GSE52583_annot.cds_Clustered.RData", compress=TRUE) #++++++++++++++++++++++++++++++
  
  ############# plot -------------
  pc <-  plot_pc_variance_explained(annot.cds, return_all = T) # norm_method='log'
  plot(pc$p)
  pc$variance_explained[1:3]*100
  #[1]  8.703065 3.243314 2.495040
  
  
  pdf(file="PCA_averageHighgene_standarded_cds_thresholded.pdf")
  plot_cell_clusters(annot.cds, color_by = 'as.factor(cell_type)')
  plot_cell_clusters(annot.cds, color_by = 'as.factor(Cluster)')
  g <- plot_cell_clusters(annot.cds, color_by = 'as.factor(cell_type)', plot=FALSE)
  g + geom_point(aes_string(color = pData(annot.cds)$cell_type, 
                            shape=pData(annot.cds)$Cluster))                             
  dev.off()
  
  
}

###################################################################
## Section 5 (OPTINAL): constructing trajectory using Monocle  ##
###################################################################
{
  annot.cds@expressionFamily@vfamily  #[1] "negbinomial.size"
  clustering_DEG_genes <- differentialGeneTest(annot.cds[which(fData(annot.cds)$use_for_ordering),],
                                               fullModelFormulaStr = '~Cluster',
                                               relative_expr = TRUE, # by default, Whether to transform expression into relative values.
                                               cores = 1)
  save(clustering_DEG_genes, file='clustering_DEG_genes.age.Cluster.RData') #!!!!!!!!!!!!
  
  ordering_genes.Cluster <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)]
  
  table(fData(annot.cds)$use_for_ordering)
  # FALSE  TRUE 
  # 12495 10359 
  annot.cds <- setOrderingFilter(annot.cds,
                                 ordering_genes = ordering_genes.Cluster)
  annot.cds <-  reduceDimension(annot.cds, method = 'DDRTree') ##### !!! there are different methods !!!
  annot.cds <- orderCells(annot.cds)
  annot.cds <- orderCells(annot.cds, root_state = GM_state(annot.cds)) 
  
 # save(annot.cds, file="GSE52583_annot.cds_Clustered_thresholded.RData", compress=TRUE) #++++++++++++++++++++++++++++++

  # ## additionally DE test for the feature selction in section 7
  # clustering_DEG_genes <- differentialGeneTest(annot.cds[x,],
  #                                              fullModelFormulaStr = '~cell_type',
  #                                              relative_expr = TRUE, # by default, Whether to transform expression into relative values.
  #                                              cores = 1)
  # save(clustering_DEG_genes, file='clustering_DEG_genes.age.cellType.RData') #!!!!!!!!!!!!!!
  
  pdf(file="pca_averageHighgene_Standardized_trajectory.pdf")
  plot_cell_trajectory(annot.cds, color_by = "age") 
  plot_cell_trajectory(annot.cds, color_by = "cell_type") 
  plot_cell_trajectory(annot.cds, color_by = "Cluster") 
  plot_cell_trajectory(annot.cds, color_by = "Pseudotime") 
  dev.off()
  
 }  

# ### backup all selected transcripts  ########
# load(file="GSE52583_annot.cds_Clustered_thresholded.RData")
# 
# dat <- exprs(annot.cds)
# dim(dat) #22854   196
# dat <- dat[which(fData(annot.cds)$use_for_ordering),]
# dim(dat)
# #[1] 10359   196
# save(dat, file="F:/projects/BioTIP/doc/2020_/Applications/input.Rdata/GSE52583/GSE52583_monocle_counts.RData")
# cli <- pData(annot.cds)
# dim(cli) # 196  21
# save(cli, file='F:/projects/BioTIP/doc/2020_/Applications/input.Rdata/GSE52583/GSE52583_cli.RData')
# df <- fData(annot.cds)[rownames(dat),]
# dim(df) # 10359    11
# save(df, file='F:/projects/BioTIP/doc/2020_/Applications/input.Rdata/GSE52583/GSE52583_GeneFeature.RData')
# rm(dat, cli, df)

######################################################
## Section 6) Global ~3k HVG section                ##
## of coding genes, lncRNAs, and miRNAs             ##
## based on variance of FPKM across all cells       ##
######################################################
{
  ## 7.1) HVG selection ---------------------
  {
    #load(file="GSE52583_annot.cds_Clustered_thresholded.RData")
    table(fData(annot.cds)$use_for_ordering)
    # FALSE  TRUE 
    # 12495 10359
    
    # select coding genes, lncRNAs, and miRNAs for the downstream analysis ------------------------
    table(fData(annot.cds)$use_for_ordering, fData(annot.cds)$gene_biotype %in% c('lncRNA', 'miRNA','protein_coding') )
    #       FALSE  TRUE
    # FALSE   380 12115
    # TRUE    108 10251
    
    ## alternatively, could select genes based on differentially expression across time or clusters
    #  x <- which(fData(annot.cds)$num_cells_expressed>= 10 & fData(annot.cds)$use_for_ordering &
    #               fData(annot.cds)$gene_biotype %in% c('lncRNA', 'miRNA','protein_coding'))
    #  length(x)  #8900
    #  
    #   
    #  load(file='clustering_DEG_genes.age.cellType.RData')
    #  selected_genes.age.cellType <-
    #    row.names(clustering_DEG_genes)[which(clustering_DEG_genes$qval<0.05)]
    #  
    #  load(file='clustering_DEG_genes.Cluster.RData')  
    #   selected_genes.Cluster <- 
    #    row.names(clustering_DEG_genes)[which(clustering_DEG_genes$qval<0.05)]
    #  HVG <- union(selected_genes.Cluster[1:4000],
    #               selected_genes.age.cellType[1:4000]
    # )
    #  if(any(is.na(HVG))) HVG <- HVG[-which(is.na(HVG))]
    #  length(HVG) # 4531
    
    x <- which(fData(annot.cds)$use_for_ordering &
                 fData(annot.cds)$gene_biotype %in% c('lncRNA', 'miRNA','protein_coding') )
    length(x) # 10251
    vars <- rowVars(log2(exprs(annot.cds[x,])+1))
    hist(vars, 100)
    table(vars>0.5)
    # FALSE  TRUE 
    # 7053  3198 
    HVG <- rownames(annot.cds)[x][which(vars>0.5)]
    
  }

  ## 7.2) export all 196 cells  -------------------------  
  {
    cli <- pData(annot.cds)
    cli$GEO.access <- rownames(cli)
    dim(cli)  # 196  22
    
    logmat <- log2(exprs(annot.cds)[HVG,]+1)
    dim(logmat)  # [1] 3198  196
    
    sce <- SingleCellExperiment(logmat)
    colData(sce) <- DataFrame(age = cli$age,
                              cellName = cli$cellName,
                              GEO.access = cli$GEO.access,
                              putative_cell_type = cli$putative_cell_type,
                              C_by_marker = cli$CellType,
                              C_Monocle.k4 = cli$Cluster,
                              Pseudotime = cli$Pseudotime,
                              Size_Factor = cli$Size_Factor)
    rowData(sce) <- fData(annot.cds)[rownames(sce),]
    names(assays(sce)) = 'logFPKM'
    dim(sce)  # 3198  196
    
    save(sce, file='BioTIP_GSE52583_robustness/sce.RData') #!!!!!!!!!!!
    
  }
  
  ## 7.3) focusing on cells along the AT2 trajectory -------------------------  
  {
    cli <- pData(annot.cds)
    cli$GEO.access <- rownames(cli)
    y <-  which(cli$age=="18.5" & cli$CellType !="AT2" )
    length(y)   # 65
    cli <- cli[-y,]; dim(cli)  # 131  21
    # remove the factor levels of 0
    tmp <- as.vector(cli$putative_cell_type)
    tmp <- factor(tmp, levels=c('BP', 'ciliated',    'Clara','AT2'       ))
    
    logmat <- log2(exprs(annot.cds)[HVG,-y]+1)
    dim(logmat)  # [1] 4531    131
    colnames(logmat) <- rownames(cli) <- cli$cellName
      
    sce <- SingleCellExperiment(logmat)
    colData(sce) <- DataFrame(age = cli$age,
                              cellName = cli$cellName,
                              GEO.access = cli$GEO.access,
                              putative_cell_type = tmp,
                              C_by_marker = cli$CellType,
                              C_Monocle.k4 = cli$Cluster,
                              Pseudotime = cli$Pseudotime,
                              Size_Factor = cli$Size_Factor)
    rowData(sce) <- fData(annot.cds)[rownames(sce),]
    names(assays(sce)) = 'logFPKM'
    dim(sce)  # 3198 131
    
    save(sce, file='BioTIP_GSE52583_robustness/AT2.sce.RData') #!!!!!!!!!!!
    
  }
}

#################################################################
## Section 7) Prepare inputs for QuanTC on 131 cells           ##
#################################################################
## cell-cell similarity matrix for 131 cells #---------------------
# refer to http://127.0.0.1:11637/library/SC3/doc/SC3.html
# needs to customize ks accordingily per dataset, the larger range the longer running time.
# In this case, ks=3:8 are tested, 
# and for the related soft-thresholding clustering (QuanTC method), 
# we had take average of the Consensus.Cluster-agreeable clustering results of k=4:10 to get cell-cell similarity matrix M
library(SC3)
files = c('AT2.sce.RData', 'sce.RData')
  for(j in files){
    load(file= paste0('BioTIP_GSE52583_robustness/',j))
    ## write into QuanTC inputs files #---------------------------
    if(j=='AT2.sce.RData') QuanTC.input.dir = "F:/projects/QuanTC/QuanTC-modified/Input/GSE52583.AT2/" else {
      QuanTC.input.dir = "F:/projects/QuanTC/QuanTC-modified/Input/GSE52583/"
    }
## 8.2) calculate cell-cell similarity and estiamte teh number of clusters ----------------------------
    {
    sce.sc3 = sce
    rowData(sce.sc3)$feature_symbol <- rownames(sce.sc3)
    # remove features with duplicated names
    any(duplicated(rowData(sce.sc3)$feature_symbol))  # F
    range(assays(sce)$logFPKM)
    # [1]     0.00000 13.97788
    
    ### to run SC3 successfully, transform sparsematrix to matrix  !!!!!!!!!!!!!
    logcounts(sce.sc3) <- as.matrix(assays(sce)$logFPKM)
    sum(logcounts(sce.sc3)<1e-16,2)/nrow(logcounts(sce.sc3))>0.95  # TRUE therefore no cell being filtered
    counts(sce.sc3) <- as.matrix(2^assays(sce)$logFPKM -1)
    
    # NOT repeat !!!! 
    # # biology: boolean parameter, defines whether to compute differentially expressed genes, marker genes and cell outliers.
    set.seed(2020)
    sce.sc3 <- sc3(sce.sc3, ks = 3:8, biology = FALSE) # svm_max = 5000 is default!!!
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
    # $ k_estimation   : num 5
    
    # to save space, transform back matrix to sparse matrix
    assayNames(sce.sc3)
    #[1] "logFPKM"   "logcounts" "counts" 
    assays(sce.sc3) <- assays(sce.sc3)[1]
    
    if(j=='AT2.sce.RData') save(sce.sc3, file='AT2.sce_SC3.RData') else save(sce.sc3, file='sce_SC3.RData') ##!!!!!!!!!!!!!!!!!!! 
    gc()
    # END DO NOT REPET !!!!!!!!!!!!!
  }   
## 8.3) writing input fiels for QuanTC -----------------------------------------
    {
    for(j in files){
      load(file= paste0('BioTIP_GSE52583_robustness/',j))
      ## write into QuanTC inputs files #---------------------------
      if(j=='AT2.sce.RData') {
        QuanTC.input.dir = "F:/projects/QuanTC/QuanTC-modified/Input/GSE52583.AT2/" 
        load(file='AT2.sce_SC3.RData')}
      else {
        QuanTC.input.dir = "F:/projects/QuanTC/QuanTC-modified/Input/GSE52583/"
        load(file='sce_SC3.RData')
      }  
  
      M_3 = (sce.sc3@metadata$sc3$consensus$`3`$consensus)
      M_5 = (sce.sc3@metadata$sc3$consensus$`5`$consensus)
      M_7 = (sce.sc3@metadata$sc3$consensus$`7`$consensus)
      # take average of the Consensus.Cluster-agreeable clustering results to get cell-cell similarity matrix M
      
      M = (M_3+M_5+M_7)/3 
      write.csv(M, file=paste0(QuanTC.input.dir,'Treutlein2014_cell-cell.csv'), row.names=F, col.names=F)   
      
      logmat <- as.matrix(assays(sce)$logFPKM)
      dim(logmat)  # [1] 4720  131
      write.table(round(logmat,4), file= paste0(QuanTC.input.dir,'Treutlein2014_log2.FPKM.txt'), row.names=FALSE, col.names=FALSE, sep='\t')
      write.table(rownames(logmat), file= paste0(QuanTC.input.dir,'Treutlein2014_gene_name.txt'), row.names=FALSE, quote=FALSE, col.names=FALSE) 
      cli <- colData(sce)
      write.table(cli$cellName %>% as.vector(), 
                  file= paste0(QuanTC.input.dir,'Treutlein2014_cell_name.txt'), row.names=FALSE, quote=FALSE, col.names=FALSE) 
      write.table(cli$C_by_marker %>% as.vector(),
                  file= paste0(QuanTC.input.dir,'Treutlein2014_CellType.txt'), row.names=FALSE, quote=FALSE, col.names=FALSE) 
      true_label <- as.vector(cli$age)
      true_label[which(true_label=='14.5')]=1   # pseudo time, numeric values
      true_label[which(true_label=='16.5')]=2  
      true_label[which(true_label=='18.5')]=3  
      true_label[which(true_label=='Adult')]=4  
      write.table(true_label, 
                  file= paste0(QuanTC.input.dir,'Treutlein2014_CellAge.txt'), row.names=FALSE, quote=FALSE, col.names=FALSE) 
      
    }
}
}
####################################################################################
## section 8) cluster 131 cells of 4531 genes using different methods             ## 
## Note that excluded cells are which(cli$age=="18.5" & cli$CellType !="AT2" )    ##
####################################################################################
setwd('F:/projects/BioTIP/result/GSE52583')

library(dplyr)
library(scater)
library(scran)
library(SC3)
library(Seurat) #Seurat version 4.0.6
#library(leiden) 
#library(monocle3)

j='AT2.sce.RData'
subDir = 'BioTIP_GSE52583_robustness/'
QuanTCDir = 'QuanTC_Output/AT2/'

{
  parameters = list()
  parameters$k = 10 # An integer number of nearest neighboring cells to use when creating the k nearest neighbor graph for Louvain/Leiden/SNNGraph clustering.
  
  ################################################################################
  ## 8.0) load the R object                                                                                         
  ################################################################################

  load(paste0(subDir,'/',j))
  sce
  # dim: 3198 131 
  # metadata(0):
  # assays(1): logFPKM
  # rownames(4531): 0610007P08Rik 0610007P22Rik ... Zyx l7Rn6
  # rowData names(11): ensembl_gene_id gene_short_name ...
  # num_cells_expressed use_for_ordering
  # colnames: NULL
  # colData names(8): age cellName ... Pseudotime Size_Factor
  # reducedDimNames(0):
  #   altExpNames(0):
    
  table(colData(sce)$C_by_marker)
  # Ambiguous       AT1       AT2   Unknown 
  #         2        13        77        39 
  table(colData(sce)$age)
  #  14.5  16.5  18.5 Adult 
  #   45    27    15    44
   
  ########################################################################################
  # 8.1) # extract SNNGraph clusters (by scran) with two settings for the parameter k
  # The parameter k indicates the number of nearest neighbors to consider during graph construction, whicc
  # we set to 5 for small number of cells (e.g., <500) and 10 to large number of sequenced cells (e.g., 4k).
  # https://nbisweden.github.io/single-cell_sib_scilifelab/session-clustering/clustering.html
  ########################################################################################
  ## Calculate size factors and normalize has been excluded becasue the downloaded FPKM has been scaled
  # sce <- scran::computeSumFactors(sce, min.mean = 0.1, assay.type ='logFPKM')
  # sce <- scater::normalize(sce)
  # logcounts(sce) <- as.matrix(logcounts(sce))
  { 
  ## Fit variance trend and apply denoising PCA
  new.trend <- scran::modelGeneVarByPoisson(x = sce, assay.type ='logFPKM')
  # means <- rowMeans(assays(sce)$logFPKM)
  # vars <- rowVars(assays(sce)$logFPKM)
  # fit <- scran::fitTrendVar(means, vars)
  # fit$trend <- new.trend
  dec <- scran::modelGeneVar(sce, assay.type ='logFPKM')
  set.seed(123)
  sce <- scran::denoisePCA(sce, technical = new.trend, assay.type ='logFPKM')
  reducedDimNames(sce)  #[1] "PCA"
   
  # k: An integer scalar specifying the number of nearest neighbors to consider during graph construction.
  SNNGraph.ID <- scran::buildSNNGraph(sce, k= parameters$k, use.dimred = 'PCA')
  SNNGraph.ID <- igraph::cluster_walktrap(SNNGraph.ID)$membership
  
  # check the agreeement between new and original clusters using the SNNGraph method
  table(as.vector(sce$age), SNNGraph.ID)
  #     SNNGraph.ID
  #        1  2  3
  # 14.5   0 40  5
  # 16.5   1  2 24
  # 18.5  12  0  3
  # Adult 44  0  0
  
  colData(sce)$C_SNNGraph_k10 = factor(SNNGraph.ID)
  
  SNNGraph.ID <- scran::buildSNNGraph(sce, k= 20, use.dimred = 'PCA') # the  same
  SNNGraph.ID <- scran::buildSNNGraph(sce, k= 8, use.dimred = 'PCA') # the  same
  SNNGraph.ID <- igraph::cluster_walktrap(SNNGraph.ID)$membership
  table(as.vector(sce$age), SNNGraph.ID)
  # SNNGraph.ID
  #        1  2  3  4  5  6  7
  # 14.5   7  0  0  0 16  0 22
  # 16.5  24  0  0  0  2  1  0
  # 18.5   2  0  0 11  0  2  0
  # Adult  0 12 22  1  0  9  0
  
  colData(sce)$C_SNNGraph_k8 = factor(SNNGraph.ID)
  
  #save(sce, file=paste0(subDir,'/',j), compress=TRUE) # !!!!!!!!!!!!!!!!!
  } 
  
  ################################################################
  # 8.2) # extract consensus clusters  
  # refer to http://127.0.0.1:11637/library/SC3/doc/SC3.html
  # needs to customize ks accordingily per dataset, the larger range the longer running time.
  # In this case, ks=3,5,7 are tested, 
  # and for the related soft-thresholding clustering (QuanTC method), 
  # we had take average of the Consensus.Cluster-agreeable clustering results of k=3,5,7 to get cell-cell similarity matrix M
  ##################################################################

  if(j=="AT2.sce.RData") load('AT2.sce_SC3.RData')  else load('sce_SC3.RData') #!!!!!!!!!!!!!!!
  {
  sce.sc3  # optimale num =5 by SC3
   
  ## manually pick the optimal matches to follow up
  table(as.vector(sce$age), colData(sce.sc3)$sc3_5_clusters)
  #       1  2  3  4  5
  # 14.5   0  0  0 37  8
  # 16.5   0  0 27  0  0
  # 18.5  14  0  0  0  1
  # Adult  0 44  0  0  0
  
  
  # load(file='sce_E8.25_HEP.RData')
  colData(sce)$C_consensus_ks3 = colData(sce.sc3)$sc3_3_clusters
  colData(sce)$C_consensus_ks5 = colData(sce.sc3)$sc3_5_clusters
  colData(sce)$C_consensus_ks7 = colData(sce.sc3)$sc3_7_clusters
  
  rm(sce.sc3)
  
  # save(sce, file=paste0(subDir,'/',j), compress=TRUE) # !!!!!!!!!!!!!!!!!
  } 
  
  ################################################################
  # 8.3) # extract Leiden clustering (using Seurat)  
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
    logcounts(sce) <- assays(sce)$logFPKM
    counts(sce) <- as.matrix(2^logcounts(sce)-1)
    
    # convert from SingleCellExperiment
    sce.seurat <- as.Seurat(sce)
    # Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
    sce.seurat
    
    # Computes the k.param nearest neighbors for a given dataset
    ndim = dim(sce.seurat[['PCA']])[2]
    sce.seurat <- FindNeighbors(sce.seurat, reduction = "PCA", k.param = parameters$k, dims = 1:ndim)
    
    sce.seurat <- FindClusters(sce.seurat, resolution = 0.4, algorithm = 4) # smaller  number of communities 
    table(Idents(sce.seurat), as.vector(colData(sce)$age))  
    #  14.5 16.5 18.5 Adult
    # 1    0    0   11    35
    # 2    6   24    2     0
    # 3   24    0    0     0
    # 4   15    2    0     0
    # 5    0    1    2     9
    colData(sce)$C_Leiden_0.4 = Idents(sce.seurat)
    
    sce.seurat <- FindClusters(sce.seurat, resolution = 0.8, algorithm = 4) # smaller  number of communities 
    table(Idents(sce.seurat), as.vector(colData(sce)$age))  
    #  14.5 16.5 18.5 Adult
    # 1    6   24    2     0
    # 2   24    0    0     0
    # 3    0    0    0    22
    # 4   15    2    0     0
    # 5    0    0   11     1
    # 6    0    0    0    12
    # 7    0    1    2     9
    colData(sce)$C_Leiden_0.8 = Idents(sce.seurat)
    
    sce.seurat <- FindClusters(sce.seurat, resolution = 1.2, algorithm = 4) # smaller  number of communities 
    table(Idents(sce.seurat), as.vector(colData(sce)$age))  
    #     14.5 16.5 18.5 Adult
    # 1   24    0    0     0
    # 2    0    0    0    22
    # 3    6   10    2     0
    # 4   15    2    0     0
    # 5    0   14    0     0
    # 6    0    0   11     1
    # 7    0    1    2     9
    # 8    0    0    0    12
    colData(sce)$C_Leiden_1.2 = Idents(sce.seurat)
    
    # In this case, remove teh pseudo counts
    counts(sce) <- logcounts(sce) <- NULL
    ## save(sce, file=paste0(subDir,'/',j), compress=TRUE) # !!!!!!!!!!!!!!!!!
    rm(sce.seurat)
    
  } 
  
  ################################################################
  # 8.4 ) # extract the QUANTC-assigned clusters
  # k=4 was optimalized by QuanTC pipeline
  # refer to QuanTC_soft.thresholding_clusters.m
  ##################################################################
  {
  C_TC <- read.table(paste0(QuanTCDir,'C_TC.txt'))
  C_TC <- C_TC[,1]
  length(C_TC) #[1] 131

  index_TC <- read.table(paste0(QuanTCDir,'index_TC.txt'))
  index_TC <- index_TC[,1]
  unique(C_TC[index_TC]) # 5 verified the C_TC is the cluster ID generated by QuanTC

  ## replace the QuanTC.cluster IDs, TC is the last, to be consistented with those shoing in Fig S2
  tmp <- data.frame(C_TC)
  tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 1, 'C1'))
  tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 2, 'C2'))
  tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 3, 'C3'))
  tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 4, 'C4'))
  tmp <- tmp %>% mutate(C_TC = replace(C_TC, C_TC == 5, 'TC'))

  table(as.vector(sce$age), tmp[,1])
  #       C1 C2 C3 C4 TC
  # 14.5  36  0  0  0  9
  # 16.5   0  0 13  0 14
  # 18.5   0  0  0  2 13
  # Adult  0 43  0  0  1

  colData(sce)$C_Soft <- tmp[,1]
}
  ## save(sce, file=paste0(subDir,'/',j), compress=TRUE) # !!!!!!!!!!!!!!!!!

  ## 8.5)  plot different clustering restuls                              
  ####################################################################
  library(scater)
  sce <- runTSNE(sce, dimred="PCA")
  {  
  x <- grep('C_', colnames(colData(sce)))
  (n=length(x)) # 11
  
  pdf(file=paste0(subDir,"/TSNE.",j,"_clustering_methods.pdf"), width=10, height=9)
  gridExtra::grid.arrange(
    plotReducedDim(sce, dimred='TSNE', colour_by='putative_cell_type', #add_legend=FALSE,
                   text_by='putative_cell_type', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('putative_cell_type')  #+ ylim(5,20)  
    ,plotReducedDim(sce, dimred='TSNE',colour_by='C_by_marker', #add_legend=FALSE,
                    text_by='C_by_marker', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('C_by gene markers') #+ ylim(5,20)
    ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Monocle.k4', #add_legend=FALSE,
                    text_by='C_Monocle.k4', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('C_Monocle.k4 with 10k genes')  #+ ylim(5,20)
    ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks3', #add_legend=FALSE,
                    text_by='C_consensus_ks3', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('C_consensus_ks3')  #+ ylim(5,20)
    ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks5', #add_legend=FALSE,
                    text_by='C_consensus_ks5', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('consensus_ks5') #+ ylim(5,20) 
    ,plotReducedDim(sce, dimred='TSNE',colour_by='C_consensus_ks7', #add_legend=FALSE,
                    text_by='C_consensus_ks7', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('consensus_ks7') #+ ylim(5,20)
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
    plotReducedDim(sce, dimred='TSNE', colour_by='C_SNNGraph_k10', #add_legend=FALSE,
                   text_by='C_SNNGraph_k10', text_size = 4, text_colour='black', point_size=0.5) +
      ggtitle('C_SNNGraph_k10') 
    ,plotReducedDim(sce, dimred='TSNE',colour_by='C_SNNGraph_k8', #add_legend=FALSE,
                    text_by='C_SNNGraph_k8', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('C_SNNGraph_k8') #+ ylim(5,20) 
    ,plotReducedDim(sce, dimred='TSNE',colour_by='C_Soft', #add_legend=FALSE,
                    text_by='C_Soft', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('soft-thresholding clusting') #+ ylim(5,20) 
    ,plotReducedDim(sce, dimred='TSNE',colour_by='age', #add_legend=FALSE,
                    text_by='age', text_size = 4, text_colour='black', point_size=0.5) + 
      ggtitle('Collection time') #+ ylim(5,20) 
    ,ncol=3, nrow=3)
  
  dev.off()
  
  }  
}

############################################################################
## scetion 9) construct and visualize the trajectory using scater package
############################################################################
fid = c('AT2.sce.RData','sce.RData')
for(k in 1:length(fid)){
  load(paste0('BioTIP_GSE52583_robustness/',fid[k]))
  
  if(!'PCA' %in% reducedDimNames(sce)){
    new.trend <- scran::modelGeneVarByPoisson(x = sce, assay.type ='logFPKM')
    dec <- scran::modelGeneVar(sce, assay.type ='logFPKM')
    set.seed(123)
    sce <- scran::denoisePCA(sce, technical = new.trend, assay.type ='logFPKM')
    reducedDimNames(sce)  #[1] "PCA"
  }
  if(!'TSNE' %in% reducedDimNames(sce)){
    sce <- runTSNE(sce, dimred="PCA")
  }
 
   ## add the label to show
  colData(sce)$label = paste(sce$age, sce$putative_cell_type, sep='_')
  
  library(scater)
  # pseudo counts
  #counts(sce) = 2^(logcounts(sce))
  
  pdf(file=paste0("BioTIP_GSE52583_robustness/trajectory_",ncol(sce),"cells.pdf"))
  
  if(fid[k]=='sce.RData') {
    by.cluster <- aggregateAcrossCells(sce, ids=sce$age, use.assay.type='logFPKM')} else {
      by.cluster <- aggregateAcrossCells(sce, ids=sce$C_Leiden_0.4, use.assay.type='logFPKM')
    }
  centroids <- reducedDim(by.cluster, "PCA")
  dmat <- dist(centroids)
  dmat <- as.matrix(dmat)
  g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
  mst <- igraph::minimum.spanning.tree(g)
  plot(mst)
  pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
  coords <- reducedDim(by.cluster, "TSNE")
  group <- rep(seq_len(nrow(pairs)), 2)
  stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)
  plotTSNE(sce, colour_by="age",
           text_by="label", text_size = 8)
  plotTSNE(sce, colour_by="age", 
           text_by="label", text_size = 8) + 
    geom_line(data=stuff, mapping=aes(x=X1, y=X2, group=group))
  
  dev.off()
  
}

k=1
load(paste0('BioTIP_GSE52583_robustness/',fid[k]))
table(sce$age, sce$C_Leiden_0.4)
#        1  2  3  4  5
# 14.5   0  6 24 15  0
# 16.5   0 24  0  2  1
# 18.5  11  2  0  0  2
# Adult 35  0  0  0  9
