
## Author:  Zhezhen Wang
## Email: zhezhen@uchicago.edu
## last update:  8/18/2019
## Acknowledgement: National Institutes of Health  R21LM012619  (Xinan H Yang)

#' Assigning Transcript Biotypes
#'
#' @aliases getBiotypes
#'
#' @description
#' The purpose of the \code{getBiotypes}() function is to class both coding and noncoding
#' transcripts into biotypes using the most recent GENCODE annotations. This
#' tool can also be used to define potential lncRNAs, given an available genome
#' transcriptome assembly (a gtf file) or any genomic loci of interest.
#'
#' @usage getBiotypes <- function(full_gr, gencode_gr, intron_gr = NULL, minoverlap = 1L)
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
#' GenomicRanges \url{https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html},
#' IRanges \url{https://www.bioconductor.org/packages/release/bioc/html/IRanges.html},
#' and BiocManager \url{http://bioconductor.org/install/index.html}.
#'
#' @return A GRanges object that returns a classified transcriptome biotypes.
#'
#' @source Refrence GRCh37 genome \url{https://www.gencodegenes.org/human/release_25lift37.html}
#' for details on gtf format visit ensemble \url{https://useast.ensembl.org/info/website/upload/gff.html}
#' @import GenomicRanges
#'
#' @references
#'
#' Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670', PMID: 30423099'
#'
#' @examples
#' # Input datasets from our package's data folder
#' library(BioTIP)
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
#' @import GenomicRanges
#' @export
#' @author Zhezhen Wang and Biniam Feleke

getBiotypes <- function(full_gr, gencode_gr, intron_gr = NULL, minoverlap = 1L) {
  require(GenomicRanges)
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
#' The \code{getReadthrough}() function is used to find long transcripts that cover more
#' than two coding regions for gene regions of interst.
#' @usage getReadthrough = function(gr,cod_gr)
#'
#' @param gr A GRanges object that shows the start and end loci on genome.
#' @param cod_gr A GRanges object contaning coding regions.
#'
#' @details For details of findOverlaps, type.partialOverlap, type.50Overlap
#'   type.toPlot, queryhits, readthrough and subjecthits see,
#'   GenomicRanges \url{https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html},
#'   IRanges \url{https://www.bioconductor.org/packages/release/bioc/html/IRanges.html},
#'    and BiocManager \url{http://bioconductor.org/install/index.html}.
#'
#' @return A GRanges object which returns overlapping regions of the classified transcript biotypes.
#'
#' @source Refrence GRCh37 genome \url{https://www.gencodegenes.org/human/release_25lift37.html}
#' for details on gtf format visit ensemble \url{https://useast.ensembl.org/info/website/upload/gff.html}
#'
#' @import GenomicRanges
#'
#' @references
#' Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic
#' score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670'
#'
#' @examples
#' #First Load datasets and libraries
#' library(BioTIP)
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
#' @export
#'
#' @author Zhezhen Wang and Biniam Feleke

getReadthrough = function(gr,cod_gr){
  require(GenomicRanges)
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
#'
#' @seealso \code{\link{sd_selection.simulation}}
#' @import psych
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

sd_selection = function(df, samplesL, cutoff = 0.01, method = 'other', control_df = NULL,control_samplesL = NULL){
  require(psych)
  if(is.null(names(samplesL))) stop('please provide name to samplesL')
  if(any(!do.call(c,lapply(samplesL,as.character)) %in% colnames(df))) stop('please check if all sample names provided in "samplesL" are in colnames of "df"')
  if(any(lengths(samplesL)<2)) stop('please make sure there are at least one sample in every state')
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
    warning('Using method "previous", make sure sampleL is in the right order')
    sds = lapply(2:ncol(test2),function(x) test2[,x]/test2[,x-1])
    tmp = tmp[-1]
    names(sds) = tmp

  }else if(method == 'itself'){
    if(cutoff>1) stop('Using method "itself", cutoff must be smaller or equal to 1')
    sds = lapply(tmp,function(x) test2[,x])
    names(sds) = tmp

  }else if(method == 'longitudinal reference'){
    if(is.null(control_df) | is.null(control_samplesL))
      stop('Using method "longitudinal reference", make sure "control_df" and "sampleL" are assigned')
    control = sapply(tmp, function(x) apply(control_df[,control_samplesL[[x]]],1,sd,na.rm = T))
    sds = lapply(tmp,function(x) test2[,x]/control[,x])
    names(sds) = tmp

  }else{
    stop("method need to be selected from 'reference','other','previous', 'itself', or 'longitudinal reference' ")
  }

  if(cutoff<1){
    topdf = nrow(df)*cutoff
    sdtop = lapply(tmp,function(x) names(sds[[x]][order(sds[[x]],decreasing = T)[1:topdf]]))
  }else{
    sdtop = lapply(tmp,function(x) names(sds[[x]][sds[[x]]>cutoff]))
  }

  names(sdtop) = tmp
  subdf = lapply(tmp,function(x) df[,as.character(samplesL[[x]])])
  names(subdf) = tmp
  subm = lapply(names(subdf), function(x) subset(subdf[[x]],row.names(subdf[[x]]) %in% sdtop[[x]]))
  names(subm)  = tmp
  return(subm)
}


#' optimization of sd selection
#'
#' @description The \code{sd_selection.simulation} filter a multi-state dataset based on a cutoff value for standard deviation per state and optimize. By default, a cutoff value of 0.01 is used. Suggested if each state contains more than 10 samples.
#'
#' @param df A dataframe of numerics, factor or character.The rows and columns
#'   represent unique transcript IDs (geneID) and sample names, respectively.
#' @param samplesL An R list of n vectors, where n equals to the number of
#'   states. Each vector gives the sample names in a state, therefore most
#'   possibly to be a vector of character or integer. Note that the vectors
#'   (sample names) has to be among the column names of the R object 'df'.
#' @param cutoff A positive numeric value. Default is 0.01, if < 1, automaticlly
#'   goes to select top x# transcripts using the a selecting method (which is
#'   either the ?reference?, ?other? stages or ?pervious?# stage), e.g. by
#'   default it will select top 1\% of the transcripts.
#' @param method Selection of methods from 'reference', 'other', 'previous',
#'   'itself', or ?longitudinal reference?. Some specific requirements for each
#'   option:
#' * 'reference', the reference has to be the first
#' * 'previous', make sure sampleL is in the right order from benign to malign
#' * 'itself', make sure the cutoff is smaller than 1
#' * 'longitudinal reference', make sure control_df and control_samplesL are not NULL
#' * default uses 'other'
#'
#' @param control_df a count matrix with unique loci row names:ID X column
#'   names:samples names of control samples, required for method 'longitudinal
#'   reference'
#' @param control_samplesL a list of characters with stages as names of control
#'   samples, required for method 'longitudinal reference'
#' @return A list of dataframe of filtered transcripts with highest standard
#'   deviation are selected from \code{df} based on a cutoff value assigned. The
#'   resulting dataframe represents a subset of the raw input \code{df}.
#' @export
#' @seealso \code{\link{sd_selection}}
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @examples
#'
#' counts = matrix(sample(1:100,30),2,30)
#' colnames(counts) = 1:30
#' row.names(counts) = paste0('loci',1:2)
#' cli = cbind(1:30,rep(c('state1','state2','state3'),each = 10))
#' colnames(cli) = c('samples','group')
#' samplesL <- split(cli[,1],f = cli[,'group'])
#' test_sd_selection <- optimize.sd_selection(counts, samplesL, B = 3, cutoff =0.01)

optimize.sd_selection = function(df,samplesL,B=100,percent=0.8,times=80,cutoff=0.01,method = 'other',control_df = NULL,control_samplesL = NULL){
   if(is.null(names(samplesL))) stop('please provide name to samplesL')
   if(any(!do.call(c,lapply(samplesL,as.character)) %in% colnames(df))) stop('please check if all sample names provided in "samplesL" are in colnames of "df"')
   if(any(lengths(samplesL)<2)) stop('please make sure there are at least one sample in every state')
   N.random = lapply(1:length(samplesL), function(x) matrix(0, nrow = nrow(df),ncol=B))
   for(i in 1:length(N.random)){
     row.names(N.random[[i]]) = row.names(df)
   }

   n = lengths(samplesL)
   k = n*percent
#   X <- nrow(counts)*top
#  Y <- sapply(lociL,nrow)

  for(i in c(1:B)) {
    random_sample = lapply(1:length(k),function(x) sample(1:n[[x]],k[[x]]))  # replace=F by default
    names(random_sample) = names(samplesL)
    selected_counts = lapply(names(samplesL), function(x) df[,random_sample[[x]]])
    test2 = sapply(selected_counts, function(x) apply(x,1,sd,na.rm = T))
    tmp = names(samplesL)
    colnames(test2) = tmp

  if(method == 'reference'){
    ref = selected_counts[[1]]
    sdref = apply(ref,1,sd,na.rm = T)
    sds = lapply(tmp,function(x) test2[,x]/sdref[,x])
    names(sds) = tmp

  }else if(method == 'other'){
    samplesL = lapply(samplesL,as.character)
    othersample = lapply(1:length(tmp), function(x) do.call(c,samplesL[-x]))
    names(othersample) = tmp
    selecteddf = do.call(cbind,selected_counts)
    #selecteds = lapply(tmp, function(x) othersample[[x]][othersample[[x]] %in% colnames(selecteddf)])
    sdother = sapply(tmp, function(x) apply(df[,othersample[[x]]],1,function(y) sd(y,na.rm = T)))

    sds = lapply(tmp,function(x) test2[,x]/sdother[,x])
    names(sds) = tmp

  }else if(method == 'previous'){
    warning('Using method "previous", make sure sampleL is in the right order')
    sds = lapply(2:ncol(test2),function(x) test2[,x]/test2[,x-1])
    names(sds) = tmp

  }else if(method == 'itself'){
    if(cutoff>1) stop('Using method "itself", cutoff must be smaller or equal to 1')
    sds = lapply(tmp,function(x) test2[,x])
    names(sds) = tmp

  }else if(method == 'longitudinal reference'){
    if(is.null(control_df) | is.null(control_samplesL))
      stop('Using method "longitudinal reference", make sure "control_df" and "sampleL" are assigned')
    control = sapply(tmp, function(x) apply(control_df[,control_samplesL[[x]]],1,sd,na.rm = T))
    sds = lapply(tmp,function(x) test2[,x]/control[,x])
    names(sds) = tmp

  }else{
    stop("method need to be selected from 'reference','other','previous','itself','longitudinal reference'")
  }

  if(cutoff<1){
    topdf = nrow(selected_counts[[1]])*cutoff
    sdtop = lapply(tmp,function(x) names(sds[[x]][order(sds[[x]],decreasing = T)[1:topdf]]))
  }else{
    sdtop = lapply(tmp,function(x) names(sds[[x]][sds[[x]]>cutoff]))
  }
  names(sdtop) = tmp
    names(N.random) = tmp
    for(j in tmp){
      N.random[[j]][sdtop[[j]],i] = 1
    }
  }

  stable = lapply(N.random,function(x) row.names(x[rowSums(x)>times,]))
  names(stable) = tmp
  subdf = lapply(tmp,function(x) df[,as.character(samplesL[[x]])])
  names(subdf) = tmp
  subm = lapply(names(subdf), function(x) subset(subdf[[x]],row.names(subdf[[x]]) %in% stable[[x]]))
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
#' @export
#' @import stringr psych
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
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

getNetwork = function(optimal,fdr = 0.05){
  require(stringr)
  require(psych)
  require(igraph)

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
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

getCluster = function(igraphL,steps = 4){
  require(igraph)
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
#' @return When method=rw: A list of \code{\link[igraph]{communities}} objects of R package
#'   igraph, whose length is the length of the input object \code{igraphL}.
#'   These \code{\link[igraph]{communities}} objects can be used for
#'   visualization when it being assigned to the 'mark.groups' parameter of the
#'   \code{\link[igraph]{plot}} function of igraph package. Otherwise this
#'   function returns a list of vectors, whose length is the length of the input
#'   object \code{igraphL}. The names of each vector are the pre-selected
#'   transcript IDs by th function \code{\link{ sd_selection}}. Each vector,
#'   whose length is the number of pre-selected transcript in a state, contains
#'   the module IDs.
#'
#' @examples
#' test = list('state1' = matrix(sample(1:10,6),3,3),'state2' =
#' matrix(sample(1:10,6),3,3),'state3' = matrix(sample(1:10,6),3,3))
#' #assign colnames and rownames to the matrix
#'
#' for(i in names(test)){
#' colnames(test[[i]]) = 1:3
#' row.names(test[[i]]) = 1:3}
#'
#' #using 'rw' or 'natural' method
#' igraphL <- getNetwork(test, fdr=1)
#' #[1] "state1:3 nodes"
#' #[1] "state2:3 nodes"
#' #[1] "state3:3 nodes"
#'
#' cl <- getCluster_methods(igraphL)
#'
#' #using 'km', 'pam' or 'hcm'
#' cl <- getCluster_methods(test, method = 'pam', cutoff=2)
#'
#' @export
#' @import igraph TSdist
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

getCluster_methods = function(igraphL, method = 'rw', cutoff = NULL){
  require(igraph)
  require(TSdist)
  require(psych)
  if(method == 'rw'){
    if(all(sapply(igraphL,class) != 'igraph'))
      stop('random walk clustering needs a list of igraph object which can be obtained using getNetwork')
    if(!is.null(cutoff)) if(cutoff%%1 !=0) warning('Please provide a integer as "cutoff" for the cluster method random walk')
    if(is.null(cutoff)) cutoff = 4
    groups = getCluster(igraphL,cutoff)
  }else if(method == 'hcm'){
    if(all(!sapply(igraphL,class) %in% c('matrix','data.frame')))
      stop('hierarchical clustering needs a list of matrix or data.frame as the 1st argument')
    if(is.null(cutoff)) stop('hierarchical clustering needs "cutoff" to be assigned as the number of clusters wanted')
    testL = lapply(igraphL, function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r)
    groupsL = lapply(1:length(testL), function(x) hclust(dist(testL[[x]]), method = "complete"))
    par(mfrow = c(1,length(groupsL)))
    sapply(groupsL, function(x) plot(x))
    groups = lapply(groupsL, function(x) cutree(x,cutoff))
  }else if(method %in% c('km','pam')){
    if(all(!sapply(igraphL,class) %in% c('matrix','data.frame')))
      stop('k-mediods or PAM clustering needs a list of matrix or data.frame as the 1st argument')
    if(is.null(cutoff)) stop('hierarchical clustering needs "cutoff" to be assigned as the number of clusters wanted')
    testL = lapply(igraphL, function(x) corr.test(t(x),adjust = 'fdr',ci=F)$r)
    groups = lapply(1:length(testL), function(x) KMedoids(testL[[x]],cutoff,distance = 'euclidean'))
  }else if(method == 'natrual'){
    warning('selecting "natural" which will not use "cutoff" parameter')
    if(all(sapply(igraphL,class) != 'igraph'))
      stop('selecting "natural" which needs a list of igraph object as the 1st argument which can be obtained using getNetwork')
    groups = lapply(1:length(igraphL), function(x) components(igraphL[[x]])$membership)

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
#' @import psych
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

getMCI = function(groups,countsL,adjust.size = FALSE){
   require(psych)
   if(all(is.na(groups))){
       warning('no loci in any of the state in the list given, please rerun getCluster_methods with a larger cutoff or provide a list of loci')
     }else{
       if(all(sapply(groups,class) =="communities")){
         membersL = lapply(groups,membership)
       }else if(any(is.na(groups)) & any(sapply(groups,class) =="communities")){
         removed = groups[is.na(groups)] ## simply remove need further modification
         groups = groups[!is.na(groups)]
         membersL = lapply(groups,membership)
       }else{
         membersL = groups
       }
       CIl = PCCol = PCCl = sdl =list()
       names(membersL) = names(groups) # probably need to be changed to names(groups) instead of names(countsL)
       if(is.null(names(groups))) warning('No names provided for "groups"')
       if(is.null(names(countsL))) warning('No names provided for "countsL"')

       loop = names(membersL)

       for(i in loop){
         test = membersL[[i]]
         if(all(is.na(test))){
           CI = sdL = PCC = PCCo = NA
         }else{
           test.counts = countsL[[i]]
           m = lapply(1:max(test),function(x) subset(test.counts, row.names(test.counts) %in% names(test[test==x])))
           comple = lapply(1:max(test),function(x) subset(test.counts, !row.names(test.counts) %in% names(test[test==x])))
           names(m) = names(comple) = 1:max(test)#letters[1:max(test)]
           #if(length(m)>26) names(m) = names(comple) = paste0(letters[1:max(test)],1:max(test))

           PCCo = lapply(names(m), function(x) abs(cor(t(comple[[x]]),t(m[[x]]))))
           PCCo_avg = sapply(PCCo,mean)

           PCC = lapply(m,function(x) abs(cor(t(x))))
           PCC_avg = sapply(PCC,function(x) (sum(x,na.rm = T)-nrow(x))/(nrow(x)^2-nrow(x)))
           sdL = lapply(m, function(x) apply(x,1,sd))
           if(adjust.size){
             CI = mapply(function(x,y,z,w) mean(x)*(y/z)*sqrt(nrow(w)), sdL,PCC_avg,PCCo_avg,m)
           }else{
             CI = mapply(function(x,y,z) mean(x)*(y/z), sdL,PCC_avg,PCCo_avg)
           }
         }
         CIl[[i]] = CI
         sdl[[i]] = sdL
         PCCl[[i]] = PCC_avg
         PCCol[[i]] = PCCo_avg
       }
       names(CIl) = names(sdl) = names(PCCl) = names(PCCol) = names(membersL)
      return(list(members = membersL,MCI = CIl,sd = sdl,PCC=PCCl,PCCo=PCCol))
    }
}

#' plot MCI barplots
#'
#' @description a barplot of MCI for all clusters in all states
#' @param MCIl a list can be obtained through getMCI
#' @param ylim a vector needed if the output barplots need to be on the same y scale
#' @param nr the number of rows to plot
#' @param nc the number of column to plot, default length(groups)
#' @param order a character vector of the order of the barplot. default NULL which is using the input list order
#' @param minsize a non-negative numeric value of minimum size allowed for a cluster
#' @param states state names should be shown on the plot, default NULL,
#' assign this if you want to show all states including states without nodes.
#' @return a list of CI score and their components
#' @export
#' @examples
#' test = list('state1' = matrix(sample(1:10,6),4,3),'state2' = matrix(sample(1:10,6),4,3),'state3' = matrix(sample(1:10,6),4,3))
#' # assign colnames and rownames to the matrix
#' for(i in names(test)){
#'  colnames(test[[i]]) = 1:3
#'  row.names(test[[i]]) = c('g1','g2','g3','g4')
#' }
#'
#' cluster = list(c(1,2,2,1),c(1,2,3,1),c(2,2,1,1))
#' names(cluster) = names(test)
#' for(i in names(cluster)){
#'  names(cluster[[i]]) = c('g1','g2','g3','g4')
#' }
#' membersL_noweight <- getMCI(cluster,test)
#' plotBar_MCI(membersL_noweight)
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

plotBar_MCI = function(MCIl,ylim = NULL,nr=1,nc = NULL,order = NULL, minsize = 3,states = NULL){
  membersL = MCIl[[1]]
  MCI = MCIl[[2]]

  if(is.null(order)){
    loop = names(membersL)
  }else{
    if(any(!order %in% names(membersL))) stop('make sure all names provided in "order" is included in names of "countsL"')
    loop = order
  }
  if(!is.null(states)) loop = states
  if(is.null(nc)) nc = length(loop)
  par(mfrow =c(nr,nc))
  for(i in loop){
    if(! i %in% names(MCI)){
     mci = m = 0
    }else{
     mci = MCI[[i]]
     m = membersL[[i]]
     nmembers = sapply(1:max(m),function(x) length(m[m == x]))
     #cex = ifelse(length(m)>20,0.7,1)
     mci = mci[!is.na(mci)]
     if(!minsize<0 & minsize != 1){
       mci = mci[!nmembers<minsize]
       nmembers = nmembers[!nmembers<minsize]
     }else{
       warning('"minisize" need to be a non')
     }
    if(length(mci) == 0) mci = 0
    }
    mci[is.na(mci)] = 0
    bar = barplot(mci,col = rainbow(length(mci), alpha = 0.3),#legend = paste0('#',names(m),'=',sapply(m,nrow)),
          main = paste0(i,' (n=', max(m),')'),ylab = 'MCI',xlab='modules',#args.legend = list(cex = cex)
          ylim =ylim,cex.axis = 1.5, cex.names = 1.5,cex.main = 1.5,cex.lab = 1.5)
    if(all(mci != 0)) text(bar,mci,nmembers,cex = 1.5)
  }
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
#' test = list('state1' = matrix(sample(1:10,6),3,3),'state2' = matrix(sample(1:10,6),3,3),'state3' = matrix(sample(1:10,6),3,3))
#' # assign colnames and rownames to the matrix
#' for(i in names(test)){
#'   colnames(test[[i]]) = 1:3
#'   row.names(test[[i]]) = 1:3}
#'
#' igraphL <- getNetwork(test, fdr=1)
#' #[1] "state1:3 nodes"
#' #[1] "state2:3 nodes"
#' #[1] "state3:3 nodes"
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
#' cl <- getCluster_methods(test,method = "pam",cutoff = 2)
#' maxCIms <- getMaxMCImember(cl, membersL_noweight[[2]], min =2)
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

getMaxMCImember = function(membersL,MCIl,minsize = 1){
  listn = names(membersL)
  if(!minsize <1){
    minsize = minsize-1
    CIl = lapply(1:length(membersL),function(x) ifelse(table(membersL[[x]])>minsize,MCIl[[x]],NA))
    module_keep = lapply(1:length(membersL), function(x) names(table(membersL[[x]])[table(membersL[[x]])>(minsize-1)]))
    membersL = lapply(1:length(membersL),function(x) membersL[[x]][membersL[[x]] %in% module_keep[[x]]])
  }else(stop('please provide a minimum size for the cluster, which should be integer that is larger than 0'))

  idx = sapply(CIl,which.max)
  maxCI =lapply(1:length(idx),function(x) names(membersL[[x]][membersL[[x]] == idx[x]]))
  names(maxCI) = listn
  names(idx) = listn
  return(list(idx = idx,members = maxCI))
}

#' Get the cluster index and network nodes of biomodule
#'
#' @description This function retrives the cluster index and network-node ids for the identified biomodule (taht showing the maximum MCI score) at each state in the study
#'
#' @param membersL A two-layer nest list of character or numeric values, any one out of the five elements outputed by the function \code{\link{getMCI}}.
#' @param idx A vector of integers which is cluster ids of biomodule (the module with the highest MCI score) per state. This is the first element of the result from \code{\link{getMaxMCImember}}
#' @return A list describing the biomodule per state, respecting to one of the five elements (members, MCI, Sd, PCC, and PCCo) by the function \code{\link{getMCI}}. The calss of the vector is dependening on the class of the input parameter \code{membersL}.
#' @export
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @seealso \code{\link{getMCI}}
#' @examples
#' test = list('state1' = matrix(sample(1:10,6),4,3),'state2' = matrix(sample(1:10,6),4,3),'state3' = matrix(sample(1:10,6),4,3))
#' # assign colnames and rownames to the matrix
#' for(i in names(test)){
#'  colnames(test[[i]]) = 1:3
#'  row.names(test[[i]]) = c('g1','g2','g3','g4')
#' }
#'
#' cluster = list(c(1,2,2,1),c(1,2,3,1),c(2,2,1,1))
#' names(cluster) = names(test)
#' for(i in names(cluster)){
#'  names(cluster[[i]]) = c('g1','g2','g3','g4')
#' }
#' membersL_noweight <- getMCI(cluster,test)
#' idx = c(1,2,1)
#' names(idx) = names(membersL_noweight[['sd']])
#' selectedSD = getMaxStats(membersL_noweight[['sd']],idx)

getMaxStats = function(membersL,idx){
  if(any(is.null(names(idx)))| any(!names(idx) %in% names(membersL))) stop('please make sure "idx" has names and all of its names is included in names of "membersL"')
  member_max = lapply(names(idx),function(x) membersL[[x]][idx[[x]]])
  names(member_max) = names(idx)
  member_max = member_max[lengths(member_max)>0]
  member_max = sapply(member_max,function(x) mean(x[[1]]))
  return(member_max)
}

## need change CI to MCI, this change is finished for functions above 6/13/2019

#' Plot the Maximized MCI per State
#'
#' @description This function generates a line plot over multiple states with the maximum MCI score per state. The module size (i.e., number of nodes) is specified at each state in parentheses.
#'
#' @param maxMCIms A list of 2 elements. The 1st element is a integer vector of module ids whose names are the state names. The 2nd element is a list of character vectors per state. The vectors are network nodes (e.g. transcript ids). This parameter can be obtained by run function \code{\link{getMaxCImember}}
#' @param MCIl A list of numeric vectors whose names are unique cluster ids. Each vector represents the MCI scores of modules in a state. This can be the second element of the output from the function \code{\link{getMCI}}.
#' @param las Numeric in {0,1,2,3}; the style of axis labels. Default is, meaning labels are parallel (=0). (link to http://127.0.0.1:21580/library/graphics/html/par.html)
#' @param order A vector of state names in the customalized order to be plotted, setting to NULL by default.
#' @param states A character vector of state names that will be shown on the plot, setting to NULL by default. Assign this if you want to show all states including a state with no resultant module. This parameter will overwrite the parameter 'order'
#' @export
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @examples
#' maxMCIms = list(c(state1 = 1,state2 = 2,state3 = 1),c(list(state1 = c('g1','g2','g3'),state2 = c('g3','g5'),state3 = c('g2','g6'))))
#' MCIl = list(state1=c('1' = 8.84,'2' = 6.4),state2 = c('1' =NA,'2' = 9.5,'3' = NA),state3 = c('1' = 2.3,'2' = 1.4))
#' plotMaxMCI(maxMCIms,MCIl)
#'

plotMaxMCI = function(maxMCIms, MCIl, las = 0, order = NULL, states = NULL){
  if(any(is.null(names(maxMCIms[[1]]))) | any(is.null(names(maxMCIms[[2]])))) stop('Please give names for the 1st and 2nd element of the "maxMCIms" as well as "MCIl"')
  if(is.null(order)){
    CI = sapply(names(maxMCIms[[1]]), function(x) MCIl[[x]][maxMCIms[[1]][[x]]])
    ln = names(maxMCIms[[1]])
  }else{
    if(any(!order %in% names(maxMCIms[[2]])))
    stop('make sure all names in "order" are in names of the 2nd element of "maxMCIms"')
    if(any(!names(maxMCIms[[2]]) %in% order))
    warning('not every state in "simulation" is plotted, make sure "order" is complete')
    CI = sapply(order, function(x) MCIl[[x]][maxMCIms[[1]][[x]]])
    ln = order
  }
  if(class(CI) == 'list'){
    warning('changing NA CI score(s) to 0')
    idx = sapply(CI,function(x) length(x)==0)
    CI[idx] = 0
    CI = do.call(c,CI)
    names(CI) = ln
  }
  if(!is.null(states)){
    CI = CI[states]
    CI[is.na(CI)] = 0
    ln = names(CI) = states
  }
  matplot(CI,type = 'l',ylab = 'MCI(m|r)',axes=F)
  len = sapply(ln,function(x) length(maxMCIms[[2]][[x]]))
  len[is.na(len)] = 0
  names(len) = ln
  text(1:length(CI),CI+0.01,paste0('(',len,')'))

  axis(2)
  axis(side=1,at=1:length(CI),labels=ln,las = las)
}

#' Obtain the identified BioTiP and its length
#' @description getCTS obtains the identified BioTiP and its length based off of MCI scores.
#' @param maxMCI A numeric vector, whose length is the number of states. This parameter is the maximum MCI score of each state, and it can be obtained from the output of \code{\link{getMaxStats}}. Names need to be included in names of \code{maxMCIms}.
#' @param maxMCIms A list of character vectors per state. The vectors are network nodes (e.g. transcript ids). This parameter is the second element of the output of the function \code{\link{getMaxMCImember}}.
#' @return A character vector, in which the elements are the unique IDs of the network nodes of the BioTiP.
#' @export
#' @examples
#' maxMCI <- c(a = 2.56, b = 8.52, c = 2.36, d = 4.81, e = 5.26)
#' maxMCIms <- list(a = c("A100", "A293", "C403"), b = c("B853", "D826", "A406"), c = c("J198", "D103", "B105"), d = c("K529", "D385", "E358"), e = c("J019", "U926", "N824"))
#' identical(names(maxMCI), names(maxMCIms))
#' # TRUE
#' getCTS(maxMCI, maxMCIms)
#' # "Length: 3"
#' # "B853" "D826" "A406"
#' @author Antonio Feliciano y Pleyto and Zhezhen Wang \email{zhezhen@@uchicago.edu}

getCTS <- function(maxMCI, maxMCIms) {
  if (is.null(names(maxMCI))) {
    stop("No names for maxMCI. Please provide names.")
  }
  if (is.null(names(maxMCIms))) {
    stop("No names for maxMCIms. Please provide names.")
  }
  if (!all(names(maxMCI) %in% names(maxMCIms))) {
    stop("Names of maxMCI has to be in maxMCIms.")
  }
  y <- maxMCIms[[names(maxMCI)[which.max(maxMCI)]]]
  print(paste0("Length: ", length(y)))
  return(y)
}

#' calculate MCI score
#'
#' @description calculate MCI score for randomly selected loci. This function will be called by simulationMCI
#'
#' @return a list MCI scores
#' @author Zhezhen Wang

getCI_inner = function(members,countsL,adjust.size){

  random_id = sample(1:nrow(countsL[[1]]),members)
  randomL = lapply(names(countsL), function(x) countsL[[x]][random_id,])
  comple = lapply(names(countsL), function(x) subset(countsL[[x]],!row.names(countsL[[x]]) %in% row.names(randomL[[x]])))
  names(randomL) = names(comple) = names(countsL)
  PCCo = lapply(names(countsL), function(x) abs(cor(t(comple[[x]]),t(randomL[[x]]))))
  PCCo_avg = sapply(PCCo,function(x) mean(x,na.rm = T))
  PCC = lapply(randomL,function(x) abs(cor(t(x))))
  PCC_avg = sapply(PCC,function(x) (sum(x,na.rm = T)-nrow(x))/(nrow(x)^2-nrow(x)))
  sdL = lapply(randomL, function(x) apply(x,1,sd))

  if(adjust.size){
    CI = mapply(function(x,y,z,w) mean(x)*(y/z)*sqrt(members), sdL,PCC_avg,PCCo_avg,members)
  }else{
    CI = mapply(function(x,y,z) mean(x)*(y/z), sdL,PCC_avg,PCCo_avg)
  }
  return(CI)
}

#' Get MCI Scores for Feature Permutation
#'
#' @description This function gets the MCI scores for randomly selected features (e.g. transcript ids)
#'
#' @param len A integer which is the length of BioTiP.
#' @param sampelL A list of vectors, whose length is the number of states. Each vector gives the sample names in a state, therefore most possibly to be a vector of character or integer. Note that the vector s (sample names) has to be among the column names of the R object 'df'.
#' @param df A numeric matrix or data framedataframe of numerics, factor or character. The rows and columns represent unique transcript IDs (geneID) and sample names, respectively
#' @param adjust.size A boolean value indicating if MCI score should be adjust by module size (the number of transcript in the module) or not. Default FALSE.
#' @param B A integer, setting the permutation with \code{B} runs. Default is 1000.
#' @param order A vector of state names in the customalized order to be plotted, setting to NULL by default.
#' @return A numeric matrix indicating the MCI scores of permutation. The dimemsion (row * column) of this matrix is the length of \code{samplesL} * \code{B}.
#' @export
#' @examples
#' counts = matrix(sample(1:100,18),3,9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1','loci2','loci3')
#' cli = cbind(1:9,rep(c('state1','state2','state3'),each = 3))
#' colnames(cli) = c('samples','group')
#' samplesL <- split(cli[,1],f = cli[,'group'])
#' simMCI = simulationMCI(2,samplesL,counts,B=2)
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

simulationMCI = function(len, samplesL, df, adjust.size = FALSE, B=1000){
  if(is.null(names(samplesL))) stop('please provide names for list countsL')
  countsL = lapply(samplesL, function(x) df[,as.character(x)])
  if(is.null(names(countsL))) names(countsL) = names(samplesL)
  m = sapply(1:B, function(x) getCI_inner(len,countsL,adjust.size))
  row.names(m) = names(countsL)
  return(m)
}

#' Simulation of Loci to Calculate the CI Score
#'
#' @description simulation of loci to calculate the Ic score
#'
#' @param MCI A named vector of max CI scores per states, can be obtained from function \code{\link{getMaxStats}}
#' @param simulation A matrix state * number of simulated times, can be obtained from function \code{\link{simulationMCI}}
#' @param las Numeric in {0,1,2,3}; the style of axis labels. Default is, meaning labels are parallel (=0). (link to http://127.0.0.1:21580/library/graphics/html/par.html)
#' @param order A vector of state names in the customalized order to be plotted, setting to NULL by default.
#' @param ylim A interger vector of length 2. Default is NULL.
#' @param main A character vector. The title of the plot. Defulalt is NULL.
#' @export
#' @examples
#' MCI = c(1:3); names(MCI) = c('a','b','c')
#' simMCI = matrix(sample(1:100,9),3,3)
#' row.names(simMCI) = names(MCI)
#' plot_MCI_Simulation(MCI,simMCI)
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

plot_MCI_Simulation = function(MCI,simulation,las = 0,order = NULL,ylim = NULL,main = NULL){
  if(is.null(names(MCI))) stop('make sure elements in "CI" have names')
  if(!is.null(order)){
    if(any(!order %in% row.names(simulation))) stop('make sure "simulation" has row.names which are in "order"')
    if(any(!row.names(simulation) %in% order)) warning('not every state in "simulation" is plotted, make sure "order" is complete')
    simulation = simulation[order,]
  }
  maxpt = max(simulation,MCI,na.rm = T)
  tmp = c(min(simulation,MCI,na.rm = T),maxpt)
  if(is.null(ylim)){
    if(min(simulation,na.rm = T)<maxpt){
      ylim = tmp
    }else{
      ylim = rev(tmp)
    }
  }
  boxplot(t(simulation),col = 'grey',ylab = 'MCI(m|r)',axes=F,ylim=ylim,main = main)

  x = which.max(MCI)
  maxCI = MCI[x]

  if(!is.null(order)){
    if(is.null(names(MCI))) stop('make sure "CI" is named using names in "order"')
  }

  axis(2)
  # customize x-axis
  if(is.null(order)){
    stages = row.names(simulation)
  }else{
    stages = order
  }
  x = which(stages == names(x))
  axis(side=1,at=1:nrow(simulation),labels=stages,las = las)
  points(x,maxCI,col = 'red',pch = 16)
}


#' get Ic score
#'
#' @description retrive Ic scores (Pearson correlation of genes / Pearson correlation of samples) for the identified BioTiP
#'
#' @param counts A  numeric matrix or data frame. The rows and columns represent unique transcript IDs (geneID) and sample names, respectively.
#' @param sampleL A list of vectors, whose length is the number of states. Each vector gives the sample names in a state, therefore most possibly to be a vector of character or integer. Note that the vector s (sample names) has to be among the column names of the R object 'df'.
#' @param genes A character vector consists of unique BioTiP ids. This can be obtained from \code{\link{getMaxMCImember}}
#' @param output A string. Please select from 'Ic','PCCg', or 'PCCs'. Uses 'Ic' by default.
#' 'PCCg' is the PCC between genes (numerator) and 'PCCs' is PCC between samples (denominator)
#' @return a list of numeric values, whose length and names are inherited from \code{sampleL}
#' @export
#' @examples
#' counts = matrix(sample(1:100,27),3,9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1','loci2','loci3')
#' cli = cbind(1:9,rep(c('state1','state2','state3'),each = 3))
#' colnames(cli) = c('samples','group')
#' samplesL <- split(cli[,1],f = cli[,'group'])
#' BioTiP = c('loci1','loci2')
#' Ic = getIc(counts,samplesL,BioTiP)
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

getIc = function(counts,sampleL,genes,output = 'Ic'){
  subsetC = subset(counts, row.names(counts) %in% genes)
  subsetC = lapply(sampleL,function(x) subsetC[,as.character(x)])

  PCCg = lapply(subsetC,function(x) abs(cor(t(x))))
  for(i in 1:length(PCCg)) PCCg[[i]][lower.tri(PCCg[[i]],diag = T)] = NA
  PCCg = sapply(PCCg,function(x) mean(x,na.rm = T))

  PCCs = lapply(subsetC,function(x) cor(x))
  for(i in 1:length(PCCs)) PCCs[[i]][lower.tri(PCCs[[i]],diag = T)] = NA
  PCCs = sapply(PCCs,function(x) mean(x,na.rm = T))

  toplot = PCCg/PCCs
  names(toplot) = names(PCCg) = names(PCCs) = names(sampleL)
  if(output == 'Ic'){
    return(toplot)
  }else if(output == 'PCCg'){
    return(PCCg)
  }else if(output == 'PCCs'){
    return(PCCs)
  }
}

#' plot a line plot of Ic score for each state
#'
#' @description plot a line plot with Ic score for each state
#'
#' @param Ic A vector with names of states. if order is not assigned then plot by the order of this vector
#' @param las Numeric in {0,1,2,3}; the style of axis labels. Default is, meaning labels are parallel (=0). (link to http://127.0.0.1:21580/library/graphics/html/par.html)
#' @param order A vector of state names in the customalized order to be plotted, setting to NULL by default.
#' @export
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @examples
#' Ic = c('state3' = 3.4,'state1' = 5.6,'state2' = 2)
#' plotIc(Ic,order = c('state1','state2','state3'))

plotIc = function(Ic,las = 0,order = NULL){
  if(!is.null(order)){
    if(any(!order %in% names(Ic))) stop('make sure "Ic" is named using names in "order"')
    if(any(!names(Ic) %in% order)) warning('not every state in "Ic" is plotted, make sure "order" is complete')
    Ic = Ic[order]
  }
  matplot(Ic,type = 'l',ylab = 'Ic',axes=F)
  axis(2)
  stages = names(Ic)
  axis(side=1,at=1:length(Ic),labels=stages,las = las)
}

#' calculating Ic scores based on simulated loci
#'
#' @description simulate \code{x} loci \code{B} to calculate the Ic score, where x should be the same as the length of identified BioTiP and B is self-defined.
#'
#' @param obs.x A integer, length of identified BioTiP
#' @param sampleL A list of vectors, whose length is the number of states. Each vector gives the sample names in a state, therefore most possibly to be a vector of character or integer. Note that the vector s (sample names) has to be among the column names of the R object 'df'.
#' @param counts A numeric matrix or dataframe in which columns are samples and rows are transcripts.
#' Each row need to have a unique row name(i.e. transcript ID)
#' @param B A integer, setting the permutation with \code{B} runs. Default is 1000.
#' @return A matrix of \code{y} rows and \code{B} columns where \code{y} is the length of sampleL and \code{B} is self-defined. Each column is a set of Ic scores calculated for each states
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @export
#' @examples
#' counts = matrix(sample(1:100,27),3,9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1','loci2','loci3')
#' cli = cbind(1:9,rep(c('state1','state2','state3'),each = 3))
#' colnames(cli) = c('samples','group')
#' samplesL <- split(cli[,1],f = cli[,'group'])
#' simulation_Ic(2,samplesL,counts,B =3)

simulation_Ic = function(obs.x,sampleL,counts,B = 1000){
  random = sapply(1:B, function(x) sample(row.names(counts),obs.x))
  tmp= sapply(1:B, function(x) getIc(counts,sampleL,random[,x],output = 'Ic'))
  row.names(tmp) = names(sampleL)
  return(tmp)
}


#' line plot of Ic score and simulated Ic scores
#'
#' @description generate a line plot of Ic score and simulated Ic scores.
#'
#' @inheritParams plotIc
#' @param simulation A numeric matrix of Ic scores in which rows are states and columns are numbers of simulated times. It can be obtained from \code{\link{simulation_loci}}
#' @export
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @examples
#' sim = matrix(sample(1:10,9),3,3)
#' row.names(sim) = paste0('state',1:3)
#' Ic = c('state1' = 3.4,'state2' = 5.6,'state3' = 2)
#' plot_Ic_Simulation(Ic,sim)

plot_Ic_Simulation = function(Ic,simulation,las = 0,ylim = NULL,order = NULL,main = NULL){
  if(any(is.null(names(Ic)))) stop('Please provide name for vector "Ic" ')

  if(!identical(names(Ic),row.names(simulation))) Ic = Ic[match(row.names(simulation),names(Ic))]
  toplot = cbind(simulation,Ic)
  if(!is.null(order)){
    if(any(!names(Ic) %in% order)) warning('not all states in Ic is plotted')
    if(any(!order %in% names(Ic))) stop('make sure "Ic" is named using names in "order"')
    toplot = toplot[order,]
  }
  matplot(toplot,type = 'l',col = c(rep('grey',ncol(toplot)-1),'red'),lty = 1,ylab = 'Ic',axes=F,ylim=ylim,main = main)
  axis(2)
  # customize x-axis
  stages = row.names(toplot)
  axis(side=1,at=1:length(stages),labels=stages,las = las)
}

#' calculating Ic scores based on simulated samples
#'
#' @description simulate \code{x} samples \code{B} times to calculate the Ic score, where x should be the same as the length of identified BioTiP and B is self-defined.
#'
#' @param counts A  numeric matrix or data frame. The rows and columns represent unique transcript IDs (geneID) and sample names, respectively.
#' @param sampleNo A integer of sample size of the tipping point
#' @param Ic A numeric value. Ic score of identified BioTiP
#' @param genes A character vector of identified BioTiP unique ids
#' @param B A integer indicating number of times to run this simulation, default 1000
#' @export
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @examples
#' counts = matrix(sample(1:100,27),3,9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1','loci2','loci3')
#' BioTiP = c('loci1','loci2')
#' plot_simulation_sample(counts,3,3.4,BioTiP,B=3)

plot_simulation_sample = function(counts,sampleNo,Ic,genes,B = 1000,
title = 'simulation of samples'){
  sampleL = lapply(1:B, function(x) sample(colnames(counts),sampleNo))
  tmp= sapply(1:B, function(x) getIc(counts,sampleL[x],genes,output = 'Ic'))
  p_v = length(tmp[tmp>Ic])/B
  den = density(tmp)
  xmin = min(Ic,den$x)
  xmax = max(Ic,den$x)
  plot(den,main = title,xlim = c(xmin,xmax))
  abline(v = Ic,col = 'red',lty = 2)
  x = max(den$x) - 0.2*diff(range(den$x))
  if(p_v == 0) p_v = paste('<',1/B)
  text(x,max(den$y),paste0('p-value: ',p_v))
  return(tmp)
}

