#'
#' Assigning Transcript Biotypes
#'
#' @aliases getBiotypes
#'
#' @description
#' The purpose of the getBiotypes function is to class both coding and noncoding
#' transcripts into biotypes using the most recent GENCODE annotations. This
#' tool can also be used to define potential lncRNAs, given an available genome
#' transcriptome assembly (a gtf file) or any genomic loci of interest.
#'
#' @param full_gr A GRanges object which contains either coding or noncoding
#'   transctipts. Each GRanges objects columns' requires a unique
#'   identifications.For further details refer to the GRanges package.
#' @param gencode_gr A GRanges object contatining a human Chr21 GENCODE reference
#'   annotations. A metadata column, "biotype", describes the transcript type.
#' @param intron_gr A GRanges object containing the coordinates of non-coding
#'   transcripts.
#' @param minoverlap An IRanges arguments which detects minimum overlap between
#'   two IRanges objects. For more information minoverlap argument refer to the
#'   IRanges package.
#'
#' @details For details of findOverlaps, type.partialOverlap, type.50Overlap
#' type.toPlot, queryhits, and subjecthits see
#' [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
#' [IRanges](https://www.bioconductor.org/packages/release/bioc/html/IRanges.html),
#' and [BiocManager](http://bioconductor.org/install/index.html).
#'
#' @return A GRanges object that returns a classified transcriptome biotypes.
#'
#' @source [Refrence GRCh37 genome](https://www.gencodegenes.org/human/release_25lift37.html)
#' for details on gtf format visit [ensemble](https://useast.ensembl.org/info/website/upload/gff.html)
#' @import GenomicRanges
#'
#' @references
#'
#' Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670', PMID: 30423099'
#'
#' @examples
#' # Input datasets from our package's data folder
#' data("gencode")
#' data("intron")
#' data("ILEF")
#'
#' # Converting datasets to GRanges object
#' gencode_gr = GRanges(gencode)
#' ILEF_gr = GRanges(ILEF)
#' cod_gr = GRanges(cod)
#' intron_gr= GRanges(intron)
#'
#' # Filtering non-coding transcripts
#' getBiotypes(ILEF_gr, gencode_gr, intron_gr)
#'
#' \dontrun{getBiotypes(intron_gr)}
#' @note
#' Replace the PATH_FILE when loading your data locally.
#'
#' @import GenomicRanges IRanges GenomeInfoDb stats
#' @importFrom stats aggregate
#' @export
#' @author Zhezhen Wang and Biniam Feleke

getBiotypes <- function(full_gr, gencode_gr, intron_gr = NULL, minoverlap = 1L) {
    if (class(full_gr) != "GRanges")
        stop("please give full_gr as a \"GRanges\" object")
    if (class(gencode_gr) != "GRanges")
        stop("pealse give gencode_gr as a \"GRanges\" object")
    if (class(intron_gr) != "GRanges" & !is.null(intron_gr))
        stop("please give intron_gr as a \"GRanges\" object")
    hits = findOverlaps(full_gr, gencode_gr, type = "within", minoverlap = minoverlap)
    full = as.data.frame(full_gr)
    full$type.fullOverlap = "de novo"
    idx = as.data.frame(mcols(full_gr[queryHits(hits)]))
    if (nrow(idx) != 0) {
        idx$biotype = as.data.frame(mcols(gencode_gr[subjectHits(hits)]))[, 1]
        idx_collapse = aggregate(as.list(idx["biotype"]), idx["Row.names"], FUN = function(X) paste(unique(X),
            collapse = ", "))
        idx_full = match(idx_collapse$Row.names, full$Row.names)
        full[idx_full, ]$type.fullOverlap = idx_collapse$biotype}
    hits = findOverlaps(full_gr, gencode_gr, minoverlap = minoverlap)
    overlaps <- pintersect(full_gr[queryHits(hits)], gencode_gr[subjectHits(hits)])
    percentOverlap <- width(overlaps)/width(gencode_gr[subjectHits(hits)])
    idx = as.data.frame(mcols(full_gr[queryHits(hits)]))
    idx$biotype = as.data.frame(mcols(gencode_gr[subjectHits(hits)]))
    idx_collapse = aggregate(as.list(idx["biotype"]), idx["Row.names"], FUN = function(X) paste(unique(X),
        collapse = ", "))
    full$type.partialOverlap = "de novo"
    idx_partial = match(idx_collapse$Row.names, full$Row.names)
    full[idx_partial, ]$type.partialOverlap = idx_collapse$biotype
    idx$percentOverlap = percentOverlap
    idx_50 = subset(idx, percentOverlap >= 0.5)
    idx_50collapse = aggregate(as.list(idx_50["biotype"]), idx_50["Row.names"], FUN = function(X) paste(unique(X),
        collapse = ", "))
    full$type.50Overlap = "de novo"
    idx_50 = match(idx_50collapse$Row.names, full$Row.names)
    full[idx_50, ]$type.50Overlap = idx_50collapse$biotype
    if (!is.null(intron_gr)) {
        hits = findOverlaps(full_gr, intron_gr)
        idx = unique(as.data.frame(mcols(full_gr[queryHits(hits)])))
        full$hasIntron = "no"
        idx_intron = match(idx$Row.names, full$Row.names)
        if (length(idx_intron) != 0)
            full[idx_intron, ]$hasIntron = "yes"
    } else (full$hasIntron = NA)
    full$type.toPlot = ifelse(full$hasIntron == "yes" & full$type.50Overlap == "protein_coding", "protein_coding_intron",
        full$type.50Overlap)
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("protein_coding", x) & grepl("antisense",
        x), "protein_coding_antisense", x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("protein_coding,", x), "protein_coding_mixed",
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl(", protein_coding", x), "protein_coding_mixed",
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("lincRNA", x), "lincRNA",
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl("antisense,", x), "antisense",
        x))
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(grepl(", antisense", x), "antisense",
        x))
    label = c("protein_coding", "protein_coding_mixed", "lincRNA", "antisense", "pseudogene, processed_pseudogene",
        "pseudogene, unprocessed_pseudogene", "de novo", "protein_coding_antisense", "protein_coding_intron",
        "miRNA")
    full$type.toPlot = sapply(full$type.toPlot, function(x) ifelse(!x %in% label, "other_noncoding",
        x))
    return(full)
}

#' Overlapping Coding Regions
#'
#' @description
#' The getReadthrough function is used to find long transcripts that cover more
#' than two coding regions for gene regions of interst.
#'
#' @param gr A GRanges object that shows the start and end loci on genome.
#' @param cod_gr A GRanges object contaning coding regions.
#'
#' @details For details of findOverlaps, type.partialOverlap, type.50Overlap
#'   type.toPlot, queryhits, readthrough and subjecthits see,
#'   [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html),
#'   [IRanges](https://www.bioconductor.org/packages/release/bioc/html/IRanges.html),
#'    and [BiocManager](http://bioconductor.org/install/index.html).
#'
#' @return A GRanges object which returns overlapping regions of the classified transcript biotypes.
#'
#' @source [Refrence GRCh37 genome](https://www.gencodegenes.org/human/release_25lift37.html)
#' for details on gtf format visit [ensemble](https://useast.ensembl.org/info/website/upload/gff.html)
#'
#' @import GenomeInfoDb GenomicRanges
#' @importFrom stats aggregate
#'
#' @references
#' Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic
#' score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670'
#'
#' @examples
#' #First Load datasets
#' data("gencode")
#' data("ILEF")
#' data("cod")
#'
#' #Asigning datasets a GRanges object
#' gencode_gr = GRanges(gencode)
#' ILEF_gr = GRanges(ILEF)
#' cod_gr = GRanges(cod)
#'
#' getReadthrough(ILEF_gr, cod_gr)
#'
#' \dontrun{getReadthrough(cod_gr)}
#'
#' @note
#' Replace the path_file when loading data locally to the data directory.
#'
#' @import GenomicRanges IRanges GenomeInfoDb
#' @importFrom stats aggregate
#' @export
#'
#' @author Zhezhen Wang and Biniam Feleke

getReadthrough = function(gr,cod_gr){
  full_table = data.frame(gr)
  overlapcount = countOverlaps(gr,cod_gr)
  completeoverlap = unique(subjectHits(findOverlaps(cod_gr,GRanges(full_table$ID),type = 'within')))
  if(length(completeoverlap) == 0){
    full_table$readthrough = ifelse(overlapcount>2,1,0)
  }else{
    full_table$readthrough = ifelse(overlapcount>2 & row.names(completeoverlap) %in% completeoverlap,1,0)
  }
  gr = GRanges(subset(full_table,readthrough ==1))
  idx = subset(full_table,readthrough==1)$ID
  overlaps = as.data.frame(findOverlaps(gr,cod_gr))
  splitoverlaps = split(overlaps,f=overlaps$queryHits)
  table(sapply(splitoverlaps,nrow)>1)
  cod_grL = sapply(splitoverlaps,function(x) cod_gr[x$subjectHits])
  overlapL = sapply(cod_grL,function(x) findOverlaps(x))
  notoverlap = sapply(overlapL,function(x) identical(queryHits(x),subjectHits(x)))
  tmp = rep(TRUE,nrow(full_table))
  tmp[full_table$readthrough==1] = notoverlap
  full_table$readthrough = ifelse(full_table$readthrough==1 & !tmp,1,0)
  return(full_table)
}
#' Selecting Highly Oscillating Transcripts
#'
#' @description \code{sd_selection} pre-selects highly oscillating transcripts
#' from the input dataset \code{df}. The dataset must contain multiple sample
#' groups (or 'states'). For each state, the function filters the dataset using
#' a cutoff value for standard deviation. The default cutoff value is 0.01
#' (i.e., higher than the top 1\% standard deviation).
#'
#' @usage sd_selection(df, samplesL, cutoff = 0.01, method = "other")
#'
#' @param df A numeric matrix or data frame. The rows and columns represent
#'   unique transcript IDs (geneID) and sample names, respectively.
#'
#' @param samplesL A list of vectors, whose length is the number of states. Each
#'   vector gives the sample names in a state, therefore most possibly to be a
#'   vector of character or integer. Note that the vector s (sample names) has
#'   to be among the column names of the R object 'df'.
#'
#' @param cutoff A positive numeric value. Default is 0.01, if < 1,
#'   automatically goes to select top x  transcripts using the a selecting
#'   method (which is either the 'reference', 'other' stages or 'previous'
#'   stage), e.g. by default it will select top 1\% of the transcripts.
#'
#' @param method Selection of methods from 'reference', 'other', 'previous',
#'   'itself' or 'longitudinal reference'.
#' * 'reference': the reference (control) state is the first (or earliest) state.
#' * 'previous': make sure sampleL is in the right order from benign to malign.
#' * 'other': all other states in the dataset are gathered as a reference control state.
#' * 'itself': make sure the cutoff is smaller than 1.
#' * 'longitudinal reference': make sure control_df and control_samplesL are not NULL.
#'
#' @return \code{sd_selection} A list of data frames, whose length is the number
#'   of states. The rows in each data frame are the filtered transcripts with
#'   highest standard deviation selected from \code{df} and based on a cutoff
#'   value assigned. Each resulting data frame represents a subset of the raw
#'   input \code{df}, with the sample ID of the same state in the column.
#'
#' @export
#'
#' @examples
#' counts = matrix(sample(1:100,18),2,9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1','loci2')
#' cli = cbind(1:9,rep(c('state1','state2','state3'),each = 3))
#' colnames(cli) = c('samples','group')
#' samplesL <- split(cli[,1],f = cli[,'group'])
#' test_sd_selection <- sd_selection(counts, samplesL, 0.01)
#' @seealso \code\link{sd_selection.simulation}
#' @author Zhezhen Wang and Biniam Feleke

sd_selection <- function(df, samplesL, cutoff = 0.01, method = 'other'){
  if(is.null(names(samplesL))) stop('please provide name to samplesL')
  tmp = names(samplesL)
  samplesL = lapply(samplesL,as.character)
  test2 = sapply(tmp, function(x) apply(df[,samplesL[[x]]],1,sd,na.rm = T))
  if(method == 'reference'){
    ref = as.character(samplesL[[1]])
    sdref = apply(df[,ref],1,sd,na.rm = T)
    sds = lapply(tmp,function(x) test2[,x]/sdref)
    names(sds) = tmp
  }else if(method == 'other'){
    othersample = lapply(1: length(samplesL), function(x) do.call(c,samplesL[-x]))
    names(othersample) = tmp
    sdother = sapply(tmp, function(x) apply(df[,as.character(othersample[[x]])],1,sd,na.rm = T))
    sds = lapply(tmp,function(x) test2[,x]/sdother[,x])
    names(sds) = tmp
  }else if(method == 'previous'){
    warning('Using method "previous", make sure samplesL is in the right order')
    sds = lapply(2:ncol(test2),function(x) test2[,x]/test2[,x-1])
    tmp = tmp[-1]
    names(sds) = tmp
  }else{
    stop("method need to be selected from 'reference','other','previous'")
  }
  if(cutoff<1){
    topdf = nrow(df)*cutoff
    sdtop = lapply(tmp,function(x) names(sds[[x]][order(sds[[x]],decreasing = T)[1:topdf]]))
  }else{
    sdtop = lapply(tmp,function(x) names(sds[[x]][sds[[x]]>cutoff]))
  }
  if(method == 'reference') tmp = tmp[-1]
  names(sdtop) = tmp
  subdf = lapply(tmp,function(x) df[,as.character(samplesL[[x]])])
  names(subdf) = tmp
  subm = lapply(names(subdf), function(x) subset(subdf[[x]],row.names(subdf[[x]]) %in% sdtop[[x]]))
  names(subm)  = tmp
  return(subm)
}

#' Building Networks of Nodes
#'
#' @description This function builds one correlation network for each state
#'   (sample group) and runs over all states. The network nodes are defined by
#'   the context of the input dataset. For transcriptomic network analysis,
#'   network nodes can be the expressed transcript IDs and network links can be
#'   the correlation coefficients. Using the Pearson Correlation Coefficiency
#'   (PCC) analysis, this function assembles a correlation network of nodes
#'   (e.g., co-expressed transcripts) for each state using the R package igraph.

#' @usage getNetwork(optimal, fdr = 0.05)

#' @param optimal An R list of x numerical data frames, where x is the number of
#'   states studied. Each data frame consists of loci with high standard
#'   deviations. This object can be obtained through \code{sd_selection}
#'   function.
#' @param fdr A numeric cutoff value for a Pearson Correlation Coefficiency
#'   (PCC) analysis. Default is 0.05. Transcripts are linked into a network if
#'   their correlations meet this PCC-significance criterion.
#'
#' @return A list of igraph objects, whose length is the length of the input
#'   object \code{optimal}. Each object is a network of correlated nodes whose
#'   PCCs met the significant criteria based on the false discovery rate (FDR)
#'   control. The length of the list is the number of states with PCC networks.
#'   If there is no PCC met the significant criteria in one state, this state
#'   will be deleted from the output.
#'
#' @examples
#' test = list('state1' = matrix(sample(1:10,6),2,3),
#'  'state2'=matrix(sample(1:10,6),2,3),
#'  'state3' = matrix(sample(1:10,6),2,3))
#'
#' for(i in names(test)){
#'   colnames(test[[i]]) = 1:3
#'   row.names(test[[i]]) = 1:2}
#'
#' igraphL <- getNetwork(test, fdr=1)
#' #[1] "state1:2 nodes"
#' #[1] "state2:2 nodes"
#' #[1] "state3:2 nodes
#'
#' @export
#' @importFrom stringr psych igraph
#' @author Zhezhen Wang and Biniam Feleke

getNetwork <- function(optimal, fdr = 0.05){
  rL = lapply(optimal,function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r)
  names(rL) = names(optimal)
  pL = lapply(optimal,function(x) corr.test(t(x),adjust = 'fdr',ci=F)$p)
  if(is.null(names(rL))) stop('give names to the input list')
  igraphL = list()
  for(i in names(rL)){
    test = rL[[i]]
    test.p = pL[[i]]
    test[lower.tri(test,diag = T)] = NA
    #test.p[lower.tri(test,diag = T)] = 1
    tmp = lapply(1:nrow(test),function(x) test[x,test.p[x,]<fdr])
    tmp_name = lapply(1:nrow(test),function(x) which(test.p[x,]<fdr))
    idx = which(lengths(tmp_name)==1)
    for(j in idx){
      names(tmp[[j]]) = names(tmp_name[[j]])
    }
    names(tmp) = row.names(test)
    edges = stack(do.call(c,tmp))
    edges = subset(edges, !is.na(values))
    tmp2 = subset(edges,grepl('\\.[1-9,A-z]\\.',ind))
    if(nrow(tmp2)!=0){
      tmp2$node1 = paste0(str_split_fixed(tmp2$ind,'\\.',3)[,1],'.',str_split_fixed(tmp2$ind,'\\.',3)[,2])
      tmp2$node2 = str_split_fixed(tmp2$ind,'\\.',3)[,3]
    }
    edges = subset(edges,!grepl('\\.[1-9,A-z]\\.',ind))
    edges$node1 = str_split_fixed(edges$ind,'\\.',2)[,1]
    edges$node2 = str_split_fixed(edges$ind,'\\.',2)[,2]
    edges = rbind(edges,tmp2)
    dim(edges) #[1] 1270    4
    dim(edges) #[1] 583   4
    edges = edges[,c('node1','node2','values')]
    edges$weight = abs(edges$values) # added in 1/8/2019
    #colnames(edges) = c('node1','node2','weight') # added in 12/18/2018
    nodes = data.frame(unique(c(edges$node1,edges$node2)))
    print(paste0(i,':',nrow(nodes),' nodes')) #[1] 48    1
    routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
    igraphL[[i]] = routes_igraph
  }
  return(igraphL)
}
#' Grouping Similar Datasets
#'
#' @description Samples with relatively similar correlation scores are grouped
#'   together. Several clusters with their respective correlation groups
#'   could be present within a networked nodes.
#'
#' @param igraphL A list of count matrix or a list of igrap objects
#' @param steps A number of steps
#' @return A list of objects with similar correlations within networks.
#'
#' @examples
#' gc <- getCluster(igraphL, steps = 4)
#' gc
#' @import igraphL
#' @author Zhezhen Wang

getCluster <- function(igraphL, steps = 4){
  #   if(is.na(step)){
  #     step = sapply(igraphL,function(x) nrow(as_data_frame(x, what="vertices"))*2)
  #     print(step)
  #   }else
  if(length(steps) == 1 & steps %% 1 == 0){               ##grepl("^[1-9]{1,}$", step) only works for 1 digit
    steps = rep(steps,length(igraphL))
  }else if(length(steps) != 1 | length(steps) != length(igraphL)){
    stop('check step: must be postive integer(s) of length 1 or length of igraphL')
  }
  groups = list()
  for(i in 1:length(igraphL)){
    if(nrow(as_data_frame(igraphL[[i]])) != 0){
      groups[[i]] = cluster_walktrap(igraphL[[i]],weight = abs(E(igraphL[[i]])$weight),steps = steps[i])
    }else{
      groups[[i]] = NA
    }
  }
  #groups = lapply(1:length(igraphL), function(x) cluster_walktrap(igraphL[[x]],weight = abs(E(igraphL[[x]])$weight),steps = steps[x])) # changed weight to abs(PCC) 12/18/2018
  names(groups) = names(igraphL)
  return(groups)
}

#' Clustering Network Nodes
#'
#' @description This function runs over all states which are grouped samples.
#'   For each state, this function splits the correlation network generated from
#'   the function \code{\link{getNetwork}} into several sub-networks (which we
#'   called 'module'). The network nodes will be defined by the end-user. For
#'   transcriptome analysis, network nodes can be the expressed transcripts. The
#'   outputs of this function include the module IDs and node IDs per module.
#'
#' @usage getCluster_methods(igraphL, method = "rw", cutoff = NULL)
#'
#' @param igraphL A list of numerical matrices or a list of igraph objects. The
#'   list of igraph objects can be the output from the getNetwork function.
#' @param steps A mathematical clustering model for analyzing network nodes.
#'   Default is a random walk ('rw'). A method could be 'rw', 'hcm', 'km',
#'   'pam', or 'natural', where:
#' * rw: random walk using cluster_walktrap function in igraph package.
#'   'igraphL' has to be a list of igraph.
#' * hcm: hierarchical clustering using function [hclust](http://127.0.0.1:17309/library/stats/html/hclust.html)
#'   and [dist](http://127.0.0.1:17309/library/stats/html/dist.html), using method
#'   'complete'.
#' * km and pam: k-medoids or PAM algorithm using [KMedoids](http://127.0.0.1:17309/library/TSdist/html/KMedoids.html).
#' * natrual: if nodes are disconnected, they may naturally cluster and form
#'   sub-networks.
#' @param cutoff A numeric value, default is NULL. For each method it means:
#' * rw: the number of steps needed, see [cluster_walktrap](https://igraph.org/r/doc/cluster_walktrap.html)
#'   for more detail. If "cutoff" is not assigned, default of 4 will be used.
#' * hcm, km and pam: number of clusters wanted. No default assigned.
#' * natrual: does not use this parameter.
#'
#' @return When method=rw: A list of \code{communities} objects of R package
#'   igraph, whose length is the length of the input object \code{igraphL}.
#'   These \code\link[communities](https://www.rdocumentation.org/packages/igraph/versions/0.4.4/topics/communities)
#'   objects can be used for visualization when it being assigned to the
#'   'mark.groups' parameter of the
#'   [plot](http://127.0.0.1:17309/library/graphics/html/plot.html) function of
#'   igraph package. Otherwise this function returns a list of vectors, whose
#'   length is the length of the input object \code{igraphL}. The names of each
#'   vector are the pre-selected transcript IDs by th function \code{\link{
#'   sd_selection}}. Each vector, whose length is the number of pre-selected
#'   transcript in a state, contains the module IDs.
#'
#' @examples
#' test = list('state1' = matrix(sample(1:10,6),2,3),'state2' =
#' matrix(sample(1:10,6),2,3),'state3' = matrix(sample(1:10,6),2,3))
#' #assign colnames and rownames to the matrix
#'
#' for(i in names(test)){
#' colnames(test[[i]]) = 1:3
#' row.names(test[[i]]) = 1:2}
#'
#' #using 'rw' or 'natural' method
#' igraphL <- getNetwork(test, fdr=1)
#' #[1] "state1:2 nodes"
#' #[1] "state2:2 nodes"
#' #[1] "state3:2 nodes"
#'
#' cl <- getCluster_methods(igraphL)
#'
#' #using 'km', 'pam' or 'hcm'
#' cl <- getCluster_methods(test, method = 'pam', cutoff=4)
#'
#' @import igraphL
#' @author Zhezhen Wang and Biniam Feleke

getCluster_methods <- function(igraphL, method = 'rw', cutoff = NULL){
  if(method == 'rw'){
    if(all(sapply(igraphL,class) != 'igraph'))
      stop('random walk clustering needs a list of igraph object which can be obtained using getNetwork')
    if(is.null(cutoff)) cutoff = 4
    if(cutoff%%1 !=0) warning('Please provide a integer as "cutoff" for the cluster method random walk')
    groups = getCluster(igraphL,cutoff)
  }else if(method == 'hcm'){
    if(all(!sapply(igraphL,class) %in% c('matrix','data.frame')))
      stop('random walk clustering needs a list of igraph object as the 1st argument which can be obtained using getNetwork')
    testL = lapply(igraphL, function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r)
    groups = lapply(1:length(testL), function(x) hclust(dist(testL[[x]]), method = "complete"))
  }else if(method %in% c('km','pam')){
    if(all(!sapply(igraphL,class) %in% c('matrix','data.frame')))
      stop('random walk clustering needs a list of igraph object as the 1st argument which can be obtained using getNetwork')
    testL = lapply(igraphL, function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r)
    groups = lapply(1:length(testL), function(x) KMedoids(testL[[x]],3,distance = 'euclidean'))
  }else if(method == 'natrual'){
    if(all(sapply(igraphL,class) != 'igraph'))
      stop('random walk clustering needs a list of igraph object as the 1st argument which can be obtained using getNetwork')
    groups = lapply(1:length(igraphL), function(x)components(igraphL[[x]])$membership)

  }else(stop('please select from "rw", "hcm","km", "pam", "natrual" as method'))
  return(groups)
}

#' Calculating MCI Scores
#'
#' @description This function calculates a module critical index (MCI) score for
#'   each module per state within a dataset. Each module is a cluster of
#'   transcripts generated from the function \code{\link{ getCluster_methods }}.
#'   Note that a dataset should contains three or more states (samples in
#'   groups).
#' @usage getMCI(groups, countsL, plot = TRUE, adjust.size = FALSE, ylim = NULL,
#'   nr = 1, nc = length(countsL), order = NULL)
#'
#' @param groups A list of elements, whose length is the member of states. The
#'   elements could be either vectors or the R object \code{communities} of the
#'   R package \code{\link{igraph}}. If a vector, it the output of the function
#'   \code{getCluster_methods}. The names of each vector are the pre-selected
#'   transcript IDs generated by the function \code{\link{ sd_selection}}. Each
#'   vector, whose length is the number of pre-selected transcript in a state,
#'   contains the module IDs. If a \code{communities} object, it can be obtained
#'   by \code{getCluster_methods} using the "rw" method. It is also an output of
#'   the function \code{\link{sd_selection}}.
#' @param countsL A list of x numeric count matrix or x data frame, where x is
#'   the member of states.
#' @param plot a boolean to decide if plot a bar plot of CI scores or not.
#'   Default TRUE.
#' @param adjust.size A boolean value indicating if MCI score should be adjusted
#'   by module size (the number of transcripts in the module) or not. Default
#'   FALSE.
#' @param ylim A vector same y scale.
#' @param nr The number of rows to plot.
#' @param nc The number of columns to plot. Default is NULL.
#' @param order The order of the barplot, By default uses the order of 'groups'.
#'
#' @return A list of five elements (members, MCI, Sd, PCC, and PCCo). Each of
#'   element is a two-layer nest list whose length is the length of the input
#'   object \code{groups}. Each internal nested list is structured according to
#'   the number of modules identified in that state.
#' * members: vectors of unique ids
#' * MCI: the MCI score
#' * sd: standard deviation
#' * PCC: Mean of pairwised Pearson Correlation Coefficient calculated among the
#' loci in a module.
#' * PCCo: Mean of pairwised Pearson Correlation Coefficient calculated between
#' the loci in a module and the loci outside that module but inside the same
#' state.
#' @export
#'
#' @examples
#' test = list('state1' = matrix(sample(1:10,6),4,3),'state2' =
#' matrix(sample(1:10,6),4,3),'state3' = matrix(sample(1:10,6),4,3))
#'
#' #assign colnames and rownames to the matrix
#' for(i in names(test)){
#'    colnames(test[[i]]) = 1:3
#'    row.names(test[[i]]) = c('g1','g2','g3','g4')}
#'
#' cluster = list(c(1,2,2,1),c(1,2,3,1),c(2,2,1,1))
#' names(cluster) = names(test)
#' for(i in names(cluster)){
#'    names(cluster[[i]]) = c('g1','g2','g3','g4')}
#'
#'membersL_noweight <- getMCI(cluster,test)
#'
#' @importFrom psych igraph
#' @author Zhezhen Wang and Biniam Feleke

getMCI <- function(groups, countsL, plot = TRUE, adjust.size = FALSE,
                   ylim = NULL, nr=1, nc = length(countsL), order = NULL){
  if(any(sapply(groups,class) == "communities")){
    membersL = lapply(groups, membership)
  }
  else{membersL = groups
  }
  CIl = PCCol = PCCl = sdl =list()
  names(membersL) = names(countsL) # probably need to be changed to names(groups) instead of names(countsL)
  if(is.null(names(countsL))) warning('No names provided for the countsL')
  if(is.null(nc)) nc = length(groups)
  if(plot) par(mfrow =c(nr,nc))

  if(is.null(order)){
    loop = 1:length(membersL)
  }else{
    loop = order
  }
  for(i in loop){
    test = membersL[[i]]
    if(all(is.na(test))){
      CI = sdL = PCC = PCCo = NA
    }
    else{
      test.counts = countsL[[i]]
      m = lapply(1:max(test),function(x) subset(test.counts, row.names(test.counts) %in% names(test[test==x])))
      comple = lapply(1:max(test),function(x) subset(test.counts, !row.names(test.counts) %in% names(test[test==x])))
      names(m) = names(comple) = letters[1:max(test)]
      #PCCo = mapply(function(x,y) abs(cor(t(x),t(y))),comple,m)
      PCCo = lapply(names(m), function(x) abs(cor(t(comple[[x]]),t(m[[x]]))))
      PCCo_avg = sapply(PCCo,mean)
      PCC = lapply(m,function(x) abs(cor(t(x))))
      #for(j in 1:length(PCC)){
      #  PCC[[j]][lower.tri(PCC[[j]],diag = T)] = NA
      #}
      PCC_avg = sapply(PCC,function(x) (sum(x,na.rm = T)-nrow(x))/(nrow(x)^2-nrow(x)))
      #PCC_avg = sapply(PCC,function(x) mean(x,na.rm = T))
      sdL = lapply(m, function(x) apply(x,1,sd))
      if(adjust.size){
        CI = mapply(function(x,y,z,w) mean(x)*(y/z)*sqrt(nrow(w)), sdL,PCC_avg,PCCo_avg,m)
      }else{
        CI = mapply(function(x,y,z) mean(x)*(y/z), sdL,PCC_avg,PCCo_avg)
      }
      if(plot) {
        tn = ifelse(is.null(order),names(membersL)[i],i)
        cex = ifelse(length(CI)>20,0.7,1) ## added 1/9/2019
        ## 3/1/2019 changed legend 'n=' to #letters=
        barplot(MCI,legend = paste0('#',names(m),'=',sapply(m,nrow)),col = rainbow(length(m), alpha = 0.3),
                main = tn,ylab = 'CI',xlab='modules',args.legend = list(cex = cex),ylim =ylim)
      }
    }
    CIl[[i]] = CI
    sdl[[i]] = sdL
    PCCl[[i]] = PCC
    PCCol[[i]] = PCCo
  }
  names(CIl) = names(PCCol) = names(PCCl) = names(sdl) = names(countsL)
  return(list(members = membersL,CI = CIl,sd = sdl,PCC=PCCl,PCCo=PCCol))
}

#' Identifying the 'Biomodule'
#'
#' @description This function reports the 'biomodule' which is the module with
#'   the maximum Module Critical Index (MCI) scores for each state. Each state
#'   can have multiple modules (groups of subnetworks derived from the function
#'   \code{\link{ getCluster_methods}}). This function runs over all states.
#'
#' @usage getMaxCImember(membersL, MCIl, minsize = 1)
#'
#' @param membersL A list of integer vectors with unique ids as names. Each
#'   vector represents the cluster number assign to that unique id. The length
#'   of this list is equal to the number of states in the study. This can be the
#'   first element of the output from function \code{getMCI} or the output from
#'   \code{getCluster_methods}, see Examples for more detail.
#' @param MCIl A list of numeric vectors with unique cluster number as names.
#'   Each vector represents the MCI scores of that module. This can be the
#'   second element of the output from function \code{getMCI}.
#' @param minsize A numerical value of the minimum module size (the number of
#'   transcript in a cluster) to output for downstream analysis.
#'
#' @return A nested list whose length is the length of the input object
#'   \code{membersL}.  Each internal list contains two objects, one object is
#'   the vector of biomodule IDs cross states, and the other object is a list of
#'   transcript IDs (each defines the biomodule per state) cross states.
#' @export
#' @examples
#' #1st option: get the input directly from getMCI function
#' test = list('state1' = matrix(sample(1:10,6),4,3),'state2' = matrix(sample(1:10,6),4,3),'state3' = matrix(sample(1:10,6),4,3))
#'
#' # assign colnames and rownames to the matrix
#' for(i in names(test)){
#'   colnames(test[[i]]) = 1:3
#'     row.names(test[[i]]) = c('g1','g2','g3','g4')}
#'
#' cluster = list(c(1,2,2,1),c(1,2,3,1),c(2,2,1,1))
#' names(cluster) = names(test)
#' for(i in names(cluster)){
#'   names(cluster[[i]]) = c('g1','g2','g3','g4')}
#'
#' membersL_noweight <- getMCI(cluster,test)
#' maxMCIms <- getMaxMCImember(membersL_noweight[[1]], membersL_noweight[[2]], min =3)
#' #The same as
#' maxMCIms <- getMaxMCImember(cluster, membersL_noweight[[2]], min =2)
#'
#' #2nd option: get the input directly from function "getCluster_methods"
#' test = list('state1' = matrix(sample(1:10,6),2,3),'state2' = matrix(sample(1:10,6),2,3),'state3' = matrix(sample(1:10,6),2,3))
#' # assign colnames and rownames to the matrix
#' for(i in names(test)){
#'   colnames(test[[i]]) = 1:3
#'   row.names(test[[i]]) = 1:2}
#'
#' igraphL <- getNetwork(test, fdr=1)
#' #[1] "state1:2 nodes"
#' #[1] "state2:2 nodes"
#' #[1] "state3:2 nodes"
#'
#' #case1: using 'rw' method by default
#' cl <- getCluster_methods(igraphL)
#' #make sure every element in list cl is a \code{communities} object
#' sapply(cl,class)
#' #       state1        state2        state3
#' #\code{communities} \code{communities} \code{communities}
#' # if not, manually remove that state and then run
#' cluster = lapply(cl, membership)
#' maxCIms <- getMaxMCImember(cluster, membersL_noweight[[2]], min =2)
#'
#' #or run function 'getMCI' and use the 1st option
#' membersL_noweight <- getMCI(cl,test)
#'
#' # case2: using methods other than the default
#' cl <- getCluster_methods(igraphL,method = "pam")
#' maxCIms <- getMaxMCImember(cl, membersL_noweight[[2]], min =2)
#'
#' @author Zhezhen Wang and Biniam Feleke

getMaxCImember <- function(membersL,MCIl,minsize = 1){
  listn = names(membersL)
  if(!minsize <1){
    minsize = minsize-1
    CIl = lapply(1:length(membersL),function(x) ifelse(table(membersL[[x]])>minsize,CIl[[x]],NA))
    # bug unsolved case     1     2     3     4     5     6
    #                      TRUE  TRUE  TRUE  FALSE  TRUE FALSE
    #CIl = lapply(1:length(membersL),function(x) CIl[[x]][table(membersL[[x]])>minsize])
    module_keep = lapply(1:length(membersL), function(x) names(table(membersL[[x]])[table(membersL[[x]])>(minsize-1)]))
    membersL = lapply(1:length(membersL),function(x) membersL[[x]][membersL[[x]] %in% module_keep[[x]]])
  }else(stop('please provide a minimum size for the cluster, which should be integer that is larger than 0'))
  idx = sapply(CIl,which.max)
  maxCI = sapply(1:length(idx),function(x) names(membersL[[x]][membersL[[x]] == idx[x]]))
  #lapply(names(memberL_weight[[1]]),function(x) names(memberL_weight[[1]][[x]][memberL_weight[[1]][[x]]==maxCIms[[1]][[x]]]))
  names(maxCI) = listn
  names(idx) = listn
  return(list(idx = idx,members = maxCI))
}


