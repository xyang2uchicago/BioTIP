
## Author:  Zhezhen Wang; Andrew Goldstein; Yuxi Sun; Xinan H Yang
## Email: zhezhen@uchicago.edu;andrewgoldstein@uchicago.edu; ysun11@uchicago.edu; xyang2@uchicago.edu
## Last update:  03/09/2022
## Acknowledgement: National Institutes of Health  R21LM012619 
## 1/15/2022 update optimize.sd_selection() to track system running progress with the packae utils, suppress the warnning in the function getMCI() by Holly Yang
## 1/04/2021 updated getIc:  allowing NA values in the countC matrix to be calculated for PCCs
## 12/02/2020 updated the following:
## Allow direct calculation of the average of PCCs which equals to an average of PCCs matrix shrunk toward 'average', 
## thus saving the computation time,
## by adding 'none' into the characteristic arguments of PCC_sample.target.
## Note that we have the two following reasons to shut off the parameter target='average' in the function cor.shrink(): 
## --Theoretically, the 'TARGET D' outlined by Schafer and Strimmer (2005) can't be fed by the estimated average. 
## --Practically, the shrinkage towards average generates an estimated matrix, 
## --whose average value remains the same as that of its observation matrix. 
## Allow a numeric argument of 'target' or 'PCC_sample.target' to shrink towards any number between [0,1].  
## Update the examples code for avg.cor.shrink() and cor.shrink() 
## Update the function cor.shrink() for the cases when Y exists. Now the same results as that only feeding a combined (X,Y). 
## 11/2020 updated the folowing: 
## Add function getNextMaxStats() 
## Revise the example code for the function getMaxStats() 
## 10/08/2020 Allowe user to chose the shrinking target for PCC_s between ('zero,'half', and 'ave') 
## 10/05/2020 Update optimize.sd_selection() 
## 9/29/2020 Updates are the following: 
## Fix bugs in the function cor.shrink() to shrink PCC_s (towards average by modifying the theory of the the Schafer-Strimmer method) 
## plot_SS_Simulation() allows the parameer xlim 
## 8/28/2020 Updates are the following: 
## Clean the bugs when calling optimize.sd_selection() with parameters method = 'previous',method = 'reference' or cuttoff = 1
## Customalize the font size of title in the function plotBar_MCI() 
## 06/19/2020 Updates are the following: 
## Add the parameter for the function getMaxMCImember() 
## Add new functions to estimate pairwise-correlation matric which is the key statistics for tipping-point identificaion. 
## Simplified the parameter using for the function getIc(); correct typo in exampling codes. 
## Correct error in the function getNetwork(). 
## Allowing getMCI to do local estimation of a correlation matrix, which runs faster than global estimation. 

#' @import utils
#' @import MASS

#' @title Assigning Transcript Biotypes
#'
#' @description
#' The purpose of the \code{getBiotypes}() function is to class both coding and noncoding 
#' transcripts into biotypes using the most recent GENCODE annotations. This 
#' tool can also be used to define potential lncRNAs, given an available genome 
#' transcriptome assembly (a gtf file) or any genomic loci of interest. 
#'
#' @param full_gr A GRanges object which contains either coding or noncoding 
#'   transcripts. Each GRanges objects' columns requires a unique 
#'   identifications. For further details refer to the GRanges package. 
#'   
#' @param gencode_gr A GRanges object containing a human Chr21 GENCODE reference 
#'   annotation. A metadata column, "biotype", describes the transcript type. 
#'   
#' @param intron_gr A GRanges object containing the coordinates of non-coding 
#'   transcripts.
#'   
#' @param minoverlap An IRanges argument which detects minimum overlap between 
#'   two IRanges objects. For more information about minoverlap argument refer to the 
#'   IRanges package.
#'
#' @details For details of findOverlaps, type.partialOverlap, type.50Overlap 
#' type.toPlot, queryhits, and subjecthits see 
#' GenomicRanges \url{https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html}, 
#' IRanges \url{https://www.bioconductor.org/packages/release/bioc/html/IRanges.html}, 
#' and BiocManager \url{http://bioconductor.org/install/index.html}.
#'
#' @return A GRanges object that returns classified transcriptome biotypes.
#'
#' @source Reference GRCh37 genome \url{https://www.gencodegenes.org/human/release_25lift37.html} 
#' for details on gtf format visit ensemble \url{https://useast.ensembl.org/info/website/upload/gff.html}
#' @importFrom GenomicRanges findOverlaps pintersect mcols width
#'
#' @references
#'
#' Wang, Z.Z., J. M. Cunningham and X. H. Yang (2018).
#' 'CisPi: a transcriptomic score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics 34(17): 664-670',  PMID: 30423099'
#'
#' @examples
#' # Input datasets from our package's data folder
#' library(GenomicRanges)
#' data("gencode")
#' data("intron")
#' data("ILEF")
#'
#' # Converting datasets to GRanges object
#' gencode_gr = GRanges(gencode)
#' ILEF_gr = GRanges(ILEF)
#' cod_gr = GRanges(cod)
#' intron_gr = GRanges(intron)
#'
#' # Filtering non-coding transcripts
#' getBiotypes(ILEF_gr,  gencode_gr,  intron_gr)
#'
#' \dontrun{getBiotypes(intron_gr)}
#' @note
#' Replace the PATH_FILE when loading your data locally.
#'
#' @importFrom GenomicRanges findOverlaps pintersect mcols width countOverlaps GRanges
#' @export
#' @author Zhezhen Wang and Biniam Feleke

getBiotypes <- function(full_gr,  gencode_gr,  intron_gr = NULL,  minoverlap = 1L) 
{
  #  require(GenomicRanges)  
  if (all(is(full_gr) != "GRanges"))
    stop("please give full_gr as a \"GRanges\" object")
  if (all(is(gencode_gr) != "GRanges"))
    stop("pealse give gencode_gr as a \"GRanges\" object")
  if (all(is(intron_gr) != "GRanges" & !is.null(intron_gr)))
    stop("please give intron_gr as a \"GRanges\" object")
  
  hits = findOverlaps(full_gr,  gencode_gr,  type = "within",  minoverlap = minoverlap)
  full = as.data.frame(full_gr)
  full$type.fullOverlap = "de novo"
  idx = as.data.frame(mcols(full_gr[queryHits(hits)]))
  if (nrow(idx) != 0) {
    idx$biotype = as.data.frame(mcols(gencode_gr[subjectHits(hits)]))[,  1]
    idx_collapse = aggregate(as.list(idx["biotype"]),  idx["Row.names"],  
                             FUN = function(X) paste(unique(X),  collapse = ",  "))
    idx_full = match(idx_collapse$Row.names,  full$Row.names)
    full[idx_full,  ]$type.fullOverlap = idx_collapse$biotype
  }
  
  hits = findOverlaps(full_gr,  gencode_gr,  minoverlap = minoverlap)
  overlaps <- pintersect(full_gr[queryHits(hits)],  gencode_gr[subjectHits(hits)])
  percentOverlap <- width(overlaps)/width(gencode_gr[subjectHits(hits)])
  idx = as.data.frame(mcols(full_gr[queryHits(hits)]))
  idx$biotype = as.data.frame(mcols(gencode_gr[subjectHits(hits)]))
  idx_collapse = aggregate(as.list(idx["biotype"]),  idx["Row.names"],  
                           FUN = function(X) paste(unique(X),  collapse = ",  "))
  full$type.partialOverlap = "de novo"
  idx_partial = match(idx_collapse$Row.names,  full$Row.names)
  full[idx_partial,  ]$type.partialOverlap = idx_collapse$biotype
  
  idx$percentOverlap = percentOverlap
  idx_50 = subset(idx,  percentOverlap >= 0.5)
  idx_50collapse = aggregate(as.list(idx_50["biotype"]),  idx_50["Row.names"],  
                             FUN = function(X) paste(unique(X),  collapse = ",  "))
  full$type.50Overlap = "de novo"
  idx_50 = match(idx_50collapse$Row.names,  full$Row.names)
  full[idx_50,  ]$type.50Overlap = idx_50collapse$biotype
  if (!is.null(intron_gr)) {
    hits = findOverlaps(full_gr,  intron_gr)
    idx = unique(as.data.frame(mcols(full_gr[queryHits(hits)])))
    full$hasIntron = "no"
    idx_intron = match(idx$Row.names,  full$Row.names)
    if (length(idx_intron) != 0)
      full[idx_intron,  ]$hasIntron = "yes"
  } else (full$hasIntron = NA)
  full$type.toPlot = ifelse(full$hasIntron ==  "yes" & full$type.50Overlap ==  "protein_coding",  
                            "protein_coding_intron", 
                            full$type.50Overlap)
  full$type.toPlot = sapply(full$type.toPlot,  
                            function(x) ifelse(grepl("protein_coding",  x) & grepl("antisense",  x), 
                                               "protein_coding_antisense",  x))
  full$type.toPlot = sapply(full$type.toPlot,  
                            function(x) ifelse(grepl("protein_coding, ",  x),  "protein_coding_mixed", 
                                               x))
  full$type.toPlot = sapply(full$type.toPlot,  
                            function(x) ifelse(grepl(",  protein_coding",  x),  "protein_coding_mixed",  x))
  full$type.toPlot = sapply(full$type.toPlot,  
                            function(x) ifelse(grepl("lincRNA",  x),  "lincRNA",  x))
  full$type.toPlot = sapply(full$type.toPlot,  
                            function(x) ifelse(grepl("antisense, ",  x),  "antisense",  x))
  full$type.toPlot = sapply(full$type.toPlot,  
                            function(x) ifelse(grepl(",  antisense",  x),  "antisense",  x))
  label = c("protein_coding",  "protein_coding_mixed",  "lincRNA",  
            "antisense",  "pseudogene,  processed_pseudogene", 
            "pseudogene,  unprocessed_pseudogene",  "de novo",  
            "protein_coding_antisense",  "protein_coding_intron", 
            "miRNA")
  full$type.toPlot = sapply(full$type.toPlot,  
                            function(x) ifelse(!x %in% label,  "other_noncoding", x))
  return(full)
}


#' @title Overlapping Coding Regions
#'
#' @description
#' The \code{getReadthrough}() function is used to find long transcripts that cover more
#' than two coding regions for gene regions of interest.
#'
#' @param gr A GRanges object that shows the start and end loci on the genome.
#' 
#' @param cod_gr A GRanges object containing coding regions.
#'
#' @details For details of findOverlaps,  type.partialOverlap, type.50Overlap
#'   type.toPlot, queryhits, readthrough and subjecthits see, 
#'   GenomicRanges \url{https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html}, 
#'   IRanges \url{https://www.bioconductor.org/packages/release/bioc/html/IRanges.html}, 
#'    and BiocManager \url{http://bioconductor.org/install/index.html}.
#'
#' @return A GRanges object that returns overlapping regions of the classified transcript biotypes.
#'
#' @source Reference GRCh37 genome \url{https://www.gencodegenes.org/human/release_25lift37.html}.
#' For details on gtf format visit ensemble \url{https://useast.ensembl.org/info/website/upload/gff.html}.
#'
#'
#' @references
#' Wang,  Z. Z.,  J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic
#' score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670'
#'
#' @examples
#' #First Load datasets and libraries
#' library(GenomicRanges)
#' data("gencode")
#' data("ILEF")
#' data("cod")
#'
#' # Assigning datasets a GRanges object
#' gencode_gr = GRanges(gencode)
#' ILEF_gr = GRanges(ILEF)
#' cod_gr = GRanges(cod)
#'
#' getReadthrough(ILEF_gr,  cod_gr)
#'
#' \dontrun{getReadthrough(cod_gr)}
#'
#' @note
#' Replace the path_file when loading data locally to the data directory.
#'
#' @export
#'
#' @author Zhezhen Wang and Biniam Feleke

getReadthrough = function(gr, cod_gr)
{
  #  require(GenomicRanges)
  full_table = data.frame(gr)
  overlapcount = countOverlaps(gr, cod_gr)
  completeoverlap = unique(subjectHits(findOverlaps(cod_gr, GRanges(full_table$ID), type = 'within')))
  if(length(completeoverlap) ==  0){
    full_table$readthrough = ifelse(overlapcount>2, 1, 0)
  }else{
    full_table$readthrough = ifelse(overlapcount>2 & row.names(completeoverlap) %in% completeoverlap, 1, 0)
  }
  gr = GRanges(subset(full_table, readthrough ==  1))
  idx = subset(full_table, readthrough ==  1)$ID
  overlaps = as.data.frame(findOverlaps(gr, cod_gr))
  splitoverlaps = split(overlaps, f = overlaps$queryHits)
  table(sapply(splitoverlaps, nrow)>1)
  cod_grL = sapply(splitoverlaps, function(x) cod_gr[x$subjectHits])
  overlapL = sapply(cod_grL, function(x) findOverlaps(x))
  notoverlap = sapply(overlapL, function(x) identical(queryHits(x), subjectHits(x)))
  tmp = rep(TRUE, nrow(full_table))
  tmp[full_table$readthrough ==  1] = notoverlap
  full_table$readthrough = ifelse(full_table$readthrough ==  1 & !tmp, 1, 0)
  return(full_table)
}

#' @title Selecting Highly Oscillating Transcripts
#'
#' @description \code{sd_selection} pre-selects highly oscillating transcripts
#' from the input dataset \code{df}. The dataset must contain multiple samples
#' groups (or 'states'). For each state,  the function filters the dataset using
#' a cutoff value for standard deviation. The default cutoff value is 0.01
#' (i.e., higher than the top 1\% standard deviation).
#'
#' @param df A numeric matrix or data frame. The rows and columns represent
#'   unique transcript IDs (geneID) and sample names, respectively.
#'
#' @param samplesL A list of vectors,  whose length is the number of states. Each
#'   vector gives the sample names in a state. Note that the vectors (sample names) must
#'   be among the column names of the R object 'df'.
#'
#' @param cutoff A positive numeric value. Default is 0.01. If <= 1, 
#'   automatically selects top x transcripts using a selecting
#'   method (which is either the \code{reference},  \code{other} stages or \code{previous}
#'   stage),  e.g. by default it will select the top 1 percentage of the transcripts.
#'
#' @param method Selection of methods from \code{reference},\code{other}, \code{previous},  
#' default uses \code{other}. Partial match enabled.
#' * \code{itself}, or \code{longitudinal reference}. Some specific requirements for each
#'   option:
#' * \code{reference},  the reference has to be the first.
#' * \code{previous},  make sure \code{sampleL} is in the right order from benign to malign.
#' * \code{itself},  make sure the cutoff is smaller than 1.
#' * \code{longitudinal reference} make sure \code{control_df} and \code{control_samplesL} are not NULL. 
#' The row number of control_df is the same as df and all transcripts in df are also in control_df.
#'
#' @param control_df A count matrix with unique loci as row names and samples names  
#' of control samples as column names,  only used for method \code{longitudinal reference}
#'   
#' @param control_samplesL A list of characters with stages as names of control
#'   samples,  required for method 'longitudinal reference'
#'
#' @return \code{sd_selection()} A list of data frames,  whose length is the number
#'   of states. The rows in each data frame are the filtered transcripts with
#'   highest standard deviation selected from \code{df} and based on an assigned cutoff
#'   value. Each resulting data frame represents a subset of the raw
#'   input \code{df},  with the sample ID of the same state in the column.
#'
#' @export
#'
#' @examples
#' counts = matrix(sample(1:100, 18), 2, 9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1', 'loci2')
#' cli = cbind(1:9, rep(c('state1', 'state2', 'state3'), each = 3))
#' colnames(cli) = c('samples', 'group')
#' samplesL <- split(cli[, 1], f = cli[, 'group'])
#' test_sd_selection <- sd_selection(counts,  samplesL,  0.01)
#'
#' @seealso \code{\link{optimize.sd_selection}}
#' @import psych
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

sd_selection = function(df,  samplesL,  cutoff = 0.01,  
                        method = c('other', 'reference', 'previous', 'itself', 'longitudinal reference'), 
                        control_df = NULL, control_samplesL = NULL)
{
  #  require(psych)
  method = match.arg(method)
  if(is.null(names(samplesL))) stop('please provide name to samplesL')
  if(any(!do.call(c, lapply(samplesL, as.character)) %in% colnames(df))) 
    stop('please check if all sample names provided in "samplesL" are in colnames of "df"')
  if(any(lengths(samplesL)<2)) stop('please make sure there are at least one sample in every state')
  
  tmp = names(samplesL)
  samplesL = lapply(samplesL, as.character)
  test2 = sapply(tmp,  function(x) apply(df[, as.character(samplesL[[x]])], 1, sd, na.rm = TRUE))
  
  if(method ==  'reference'){
    ref = as.character(samplesL[[1]])
    sdref = apply(df[, ref], 1, sd, na.rm = TRUE)
    #sds = lapply(tmp, function(x) test2[, x]/sdref)
    sds = sapply(tmp, function(x) test2[, x]/sdref)  # corrected on 8/28/2020
    names(sds) = tmp
    
  }else if(method ==  'other'){
    othersample = lapply(1: length(samplesL),  function(x) do.call(c, samplesL[-x]))
    names(othersample) = tmp
    sdother = sapply(tmp,  
                     function(x) apply(df[, as.character(othersample[[x]])], 1, sd, na.rm = TRUE))
    
    sds = lapply(tmp, function(x) test2[, x]/sdother[, x])
    names(sds) = tmp
    
  }else if(method ==  'previous'){
    warning('Using method "previous",  make sure sampleL is in the right order')
    sds = lapply(2:ncol(test2), function(x) test2[, x]/test2[, x-1])
    tmp = tmp[-1]
    names(sds) = tmp
    
  }else if(method ==  'itself'){
    if(cutoff>1) stop('Using method "itself",  cutoff must be smaller or equal to 1')
    sds = lapply(tmp, function(x) test2[, x])
    names(sds) = tmp
    
  }else if(method ==  'longitudinal reference'){
    if(is.null(control_df) | is.null(control_samplesL))
      stop('Using method "longitudinal reference",  
           make sure "control_df" and "sampleL" are assigned')
    if(nrow(df) != nrow(control_df) | !all(row.names(df) %in% row.names(control_df)))
      stop('please make sure the row numbers of "control_df" 
           is the same as "df" and all transcripts in "df" are also in "control_df".')
    control = sapply(tmp,  
                     function(x) apply(control_df[, as.character(control_samplesL[[x]])], 
                                       1,  sd,  na.rm = TRUE))
    sds = lapply(tmp, function(x) test2[, x]/control[, x])
    names(sds) = tmp
    
  }else{
    stop("method need to be selected from 'reference', 'other', 'previous',  'itself',  
         or 'longitudinal reference' ")
  }
  
  if(cutoff<= 1){  # add = 8/28/2020
    topdf = nrow(df)*cutoff
    sdtop = lapply(tmp, function(x) names(sds[[x]][order(sds[[x]], decreasing = TRUE)[1:topdf]]))
  }else{
    sdtop = lapply(tmp, function(x) names(sds[[x]][sds[[x]]>cutoff]))
  }
  
  names(sdtop) = tmp
  subdf = lapply(tmp, function(x) df[, as.character(samplesL[[x]])])
  names(subdf) = tmp
  subm = lapply(names(subdf),  function(x) subset(subdf[[x]], row.names(subdf[[x]]) %in% sdtop[[x]]))
  names(subm) = tmp
  
  if(any(is.na(subm))) {  ## added on 2/28/2020
    a <- apply(subm[[i]], 1, function(x) sum(x,na.rm = TRUE))
    tmp <- which(is.na(a))
    if(length(tmp)>0) subm[[i]] <- subm[[i]][-a,]
    b <- apply(subm[[i]], 2, function(x) sum(x,na.rm = TRUE))
    tmp <- which(is.na(b))
    if(length(tmp)>0) subm[[i]] <- subm[[i]][,-b]
  }
  
  return(subm)
}


#' @title Optimization of sd selection
#'
#' @description The \code{optimize.sd_selection} filters a multi-state dataset 
#' based on a cutoff value for standard deviation per state and optimizes. 
#' By default, a cutoff value of 0.01 is used. Suggested if each state contains more than 10 samples.
#'
#' @param df A dataframe of numerics. The rows and columns
#'   represent unique transcript IDs (geneID) and sample names, respectively.
#'   
#' @param samplesL A list of n vectors,  where n is the number of
#'   states. Each vector gives the sample names in a state. Note that the vectors
#'   (sample names) has to be among the column names of the R object 'df'.
#'   
#' @param B An integer indicating the number of times to run this optimization, default 1000.
#' @param percent A numeric value indicating the percentage of samples will 
#' be selected in each round of simulation.
#' 
#' @param times A numeric value indicating the percentage of \code{B} times a transcript 
#' need to be selected in order to be considered a stable signature.
#' 
#' @param cutoff A positive numeric value. Default is 0.01. If <= 1, automatically
#'   goes to select top x# transcripts using a selecting method (which is
#'   either the \code{reference}, \code{other} or \code{previous} stage), e.g. by
#'   default it will select the top 1 percentage of the transcripts.
#'   
#' @param method Selection of methods from \code{reference}, \code{other}, \code{previous},  
#' default uses \code{other}. Partial match enabled.
#' * \code{itself}, or \code{longitudinal reference}. Some specific requirements for each
#'   option:
#' * \code{reference}, the reference has to be the first.
#' * \code{previous}, make sure \code{sampleL} is in the right order from benign to malign.
#' * \code{itself}, make sure the cutoff is smaller than 1.
#' * \code{longitudinal reference} make sure control_df and control_samplesL are not NULL. 
#' The row numbers of control_df is the same as df and all transcript in df are also in control_df.
#'
#' @param control_df A count matrix with unique loci as row names and samples names 
#' of control samples as column names,  only used for method \code{longitudinal reference}.
#' 
#' @param control_samplesL A list of characters with stages as names of control
#'   samples,  required for method 'longitudinal reference'.
#'   
#' @return A list of dataframe of filtered transcripts with the highest standard
#'   deviation are selected from \code{df} based on a cutoff value assigned. The
#'   resulting dataframe represents a subset of the raw input \code{df}.
#'   
#' @export
#' @seealso \code{\link{sd_selection}}
#' 
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @examples
#'
#' counts = matrix(sample(1:100, 30), 2, 30)
#' colnames(counts) = 1:30
#' row.names(counts) = paste0('loci', 1:2)
#' cli = cbind(1:30, rep(c('state1', 'state2', 'state3'), each = 10))
#' colnames(cli) = c('samples', 'group')
#' samplesL <- split(cli[, 1], f = cli[, 'group'])
#' test_sd_selection <- optimize.sd_selection(counts,  samplesL,  B = 3,  cutoff = 0.01)

optimize.sd_selection = function(df,  samplesL,  B = 100,  percent = 0.8,  
                                 times = 0.8, cutoff = 0.01, 
                                 method = c('other', 'reference', 'previous', 'itself', 'longitudinal reference'), 
                                 control_df = NULL, control_samplesL = NULL)
{
  require(utils) 
   
  method = match.arg(method)
  if(is.null(names(samplesL))) stop('please provide name to samplesL')
  if(any(!do.call(c, lapply(samplesL, as.character)) %in% colnames(df))) 
    stop('please check if all sample names provided in "samplesL" are in colnames of "df"')
  if(any(lengths(samplesL)<2)) stop('please make sure there are at least one sample in every state')
  N.random = lapply(seq_along(samplesL),  function(x) matrix(0,  nrow = nrow(df), ncol = B))
  for(i in seq_along(N.random)){
    row.names(N.random[[i]]) = row.names(df)
  }
  
  N = lengths(samplesL)
  k = round(N*percent)
  #   X <- nrow(counts)*top
  #  Y <- sapply(lociL, nrow)
  
  for(i in c(1:B)) {
    #Sys.sleep(0.5); 
    setTxtProgressBar(pb, i)
    
    random_sample_id = lapply(seq_along(k), 
                              function(x) sample(1:N[[x]], k[[x]]))  # replace = FALSE by default
    names(random_sample_id) = names(samplesL)
    selected_matrix = lapply(names(samplesL),  
                             function(x) df[, samplesL[[x]][random_sample_id[[x]]]])  ## update 10/05/2020
    test2 = sapply(selected_matrix,  
                   function(x) apply(x, 1, sd, na.rm = TRUE)) # a matrix of gene_sd in row and class in column
    tmp = names(samplesL)
    colnames(test2) = tmp
    
    if(method ==  'reference'){
      #ref = selected_counts[[1]]
      #sdref = apply(ref, 1, sd, na.rm = TRUE)
      sdref = test2[,1]  ## simplified on 8/28/2020
      #sds = lapply(tmp, function(x) test2[, x]/sdref[, x])
      sds = sapply(tmp, function(x) test2[, x]/sdref) ## correct on 8/28/2020
      names(sds) = tmp
      
    }else if(method ==  'other'){
      samplesL = lapply(samplesL, as.character)
      othersample = lapply(seq_along(tmp),  
                           function(x) do.call(c, samplesL[-x]))
      names(othersample) = tmp
      selecteddf = do.call(cbind, selected_matrix)
      #selecteds = lapply(tmp,  function(x) othersample[[x]][othersample[[x]] %in% colnames(selecteddf)])
      sdother = sapply(tmp,  
                       function(x) apply(df[, othersample[[x]]],  1,  
                                         function(y) sd(y, na.rm = TRUE)))
      
      sds = lapply(tmp, function(x) test2[, x]/sdother[, x])
      names(sds) = tmp
      
    }else if(method ==  'previous'){
      cat('Using method "previous",  make sure sampleL is in the right order')
      sds = lapply(2:ncol(test2), function(x) test2[, x]/test2[, x-1])
      tmp <- tmp[-1]  ## corrected 08/28/2020 
      names(sds) = tmp  
      
    }else if(method ==  'itself'){
      if(cutoff>1) stop('Using method "itself",  cutoff must be smaller or equal to 1')
      sds = lapply(tmp, function(x) test2[, x])
      names(sds) = tmp
      
    }else if(method ==  'longitudinal reference'){
      if(is.null(control_df) | is.null(control_samplesL))
        stop('Using method "longitudinal reference",  
             make sure "control_df" and "sampleL" are assigned')
      if(nrow(df) != nrow(control_df) | !all(row.names(df) %in% row.names(control_df)))
        stop('please make sure the row numbers of "control_df" 
             is the same as "df" and all transcripts in "df" are also in "control_df".')
      control = sapply(tmp,  function(x) 
        apply(control_df[, as.character(control_samplesL[[x]])], 1, sd, na.rm = TRUE))
      sds = lapply(tmp, function(x) test2[, x]/control[, x])
      names(sds) = tmp
      
    }else{
      stop("method need to be selected from 'reference', 
           'other', 'previous', 'itself', 'longitudinal reference'")
    }
    
    if(cutoff<= 1){
      topdf = nrow(selected_matrix[[1]])*cutoff  # each selected_counts[[i]] has the same rownames as the input df
      ## topdf = nrow(selected_matrix[[i]])*cutoff  # error-causing when (i in 1:B) larger than the lengths of samplesL
      sdtop = lapply(tmp,  
                     function(x) names(sds[[x]][order(sds[[x]], decreasing = TRUE)[1:topdf]]))
    }else{
      sdtop = lapply(tmp,  function(x) names(sds[[x]][sds[[x]]>cutoff]))
    }
    
    # cat(i, '\t')
    names(sdtop) = tmp
    names(N.random) = tmp
    for(j in tmp){
      N.random[[j]][sdtop[[j]], i] = 1 # mark the selection
    }
  }
  Sys.sleep(0.01)
  close(pb)
                     
  times = times*B
  stable = lapply(N.random, function(x) row.names(x[rowSums(x)>times, ]))
  names(stable) = tmp
  subdf = lapply(tmp, function(x) df[, as.character(samplesL[[x]])])
  names(subdf) = tmp
  subm = lapply(names(subdf),  
                function(x) subset(subdf[[x]], row.names(subdf[[x]]) %in% stable[[x]]))
  names(subm) = tmp
  
  return(subm)
}


getCluster = function(igraphL, steps = 4)
{
  #  require(igraph)
  if(length(steps) ==  1 & steps %% 1 ==  0){  ##grepl("^[1-9]{1, }$",  step) only works for 1 digit
    steps = rep(steps, length(igraphL))
  }else if(length(steps) != 1 | length(steps) != length(igraphL)){
    stop('check step: must be postive integer(s) of length 1 or length of igraphL')
  }
  groups = list()
  for(i in seq_along(igraphL)){
    if(nrow(as_data_frame(igraphL[[i]])) != 0){
      groups[[i]] = cluster_walktrap(igraphL[[i]], 
                                     weight = abs(E(igraphL[[i]])$weight), 
                                     steps = steps[i])
    }else{
      groups[[i]] = NA
    }
  }
  #groups = lapply(seq_along(igraphL),  
  #      function(x) cluster_walktrap(igraphL[[x]], weight = abs(E(igraphL[[x]])$weight), steps = steps[x])) 
  # changed weight to abs(PCC) 12/18/2018
  
  names(groups) = names(igraphL)
  return(groups)
}

#' @title Clustering Network Nodes
#'
#' @description This function runs over all states which are grouped samples.
#'   For each state,  this function splits the correlation network generated from
#'   the function \code{\link{getNetwork}} into several sub-networks (which we
#'   called 'module'). The network nodes will be defined by the end-user. For
#'   transcriptome analysis,  network nodes can be the expressed transcripts. The
#'   outputs of this function include the module IDs and node IDs per module.
#'
#' @param igraphL A list of numerical matrices or a list of igraph objects. The
#'   list of igraph objects can be the output from the getNetwork function.
#'   
#' @param method A mathematical clustering model for analyzing network nodes.
#'   Default is a random walk ('rw'). A method could be 'rw',  'hcm',  'km', 
#'   'pam',  or 'natural',  where:
#' * rw: random walk using cluster_walktrap function in igraph package.
#'   'igraphL' has to be a list of igraph.
#' * hcm: hierarchical clustering using function \link[stats]{hclust})
#'   and \link[stats]{dist},  using method
#'   'complete'.
#' * km and pam: k-medoids or PAM algorithm using \link[TSdist]{KMedoids}.
#' * natural: if nodes are disconnected,  they may naturally cluster and form
#'   sub-networks.
#'   
#' @param cutoff A numeric value,  default is NULL. For each method it means:
#' * rw: the number of steps needed, see \link[igraph]{cluster_walktrap}
#'   for more detail. If "cutoff" is not assigned,  default of 4 will be used.
#' * hcm,  km and pam: number of clusters wanted. No default assigned.
#' * natural: does not use this parameter.
#'
#' @return When method = rw: A list of \code{\link[igraph]{communities}} objects of R package
#'   igraph, whose length is the length of the input object \code{igraphL}.
#'   These \code{\link[igraph]{communities}} objects can be used for
#'   visualization when being assigned to the 'mark.groups' parameter of the
#'  \code{\link[igraph]{plot.igraph}} function of the igraph package. Otherwise this
#'   function returns a list of vectors, whose length is the length of the input
#'   object \code{igraphL}. The names of each vector are the pre-selected
#'   transcript IDs by th function \code{\link{sd_selection}}. Each vector, 
#'   whose length is the number of pre-selected transcript in a state,  contains
#'   the module IDs.
#'
#' @examples
#' test = list('state1' = matrix(sample(1:10, 6), 3, 3), 'state2' = 
#' matrix(sample(1:10, 6), 3, 3), 'state3' = matrix(sample(1:10, 6), 3, 3))
#' #assign colnames and rownames to the matrix
#'
#' for(i in names(test)){
#' colnames(test[[i]]) = 1:3
#' row.names(test[[i]]) = 1:3}
#'
#' #using 'rw' or 'natural' method
#' igraphL <- getNetwork(test,  fdr = 1)
#' #[1] "state1:3 nodes"
#' #[1] "state2:3 nodes"
#' #[1] "state3:3 nodes"
#'
#' cl <- getCluster_methods(igraphL)
#'
#' #using 'km',  'pam' or 'hcm'
#' cl <- getCluster_methods(test,  method = 'pam',  cutoff = 2)
#'
#' @export
#' @import igraph cluster TSdist
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

getCluster_methods = function(igraphL,  method = c('rw', 'hcm', 'km', 'pam', 'natural'),  cutoff = NULL)
{
  #  require(igraph)
  #  require(TSdist)
  #  require(psych)
  method <- match.arg(method)
  if(method ==  'rw'){
    if(all(sapply(igraphL, class) != 'igraph'))
      stop('random walk clustering needs a list of igraph object 
           which can be obtained using getNetwork')
    if(!is.null(cutoff)) if(cutoff%%1 != 0) 
      warning('Please provide a integer as "cutoff" for the cluster method random walk')
    if(is.null(cutoff)) cutoff = 4
    groups = getCluster(igraphL, cutoff)
  } else if(method ==  'hcm'){
    if(all(!sapply(igraphL, class) %in% c('matrix', 'data.frame')))
      stop('hierarchical clustering needs a list of matrix or data.frame as the 1st argument')
    if(is.null(cutoff)) stop('hierarchical clustering needs "cutoff" 
                             to be assigned as the number of clusters wanted')
    testL = lapply(igraphL,  function(x) corr.test(t(x), adjust = 'fdr', ci = FALSE)$r)
    groupsL = lapply(seq_along(testL),  function(x) hclust(dist(testL[[x]]),  method = "complete"))
    par(mfrow = c(1, length(groupsL)))
    sapply(groupsL,  function(x) plot(x))
    groups = lapply(groupsL,  function(x) cutree(x, cutoff))
  } else if(method %in% c('km', 'pam')){
    if(all(!sapply(igraphL, class) %in% c('matrix', 'data.frame')))
      stop('k-mediods or PAM clustering needs a list of matrix or data.frame as the 1st argument')
    if(is.null(cutoff)) stop('hierarchical clustering needs "cutoff" 
                             to be assigned as the number of clusters wanted')
    testL = lapply(igraphL,  function(x) corr.test(t(x), adjust = 'fdr', ci = FALSE)$r)
    groups = lapply(seq_along(testL),  
                    function(x) pam(testL[[x]], cutoff, metric = 'euclidean')$clustering)
  }else if(method ==  'natrual'){
    warning('selecting "natural" which will not use "cutoff" parameter')
    if(all(sapply(igraphL, class) != 'igraph'))
      stop('selecting "natural" which needs a list of igraph object 
           as the 1st argument which can be obtained using getNetwork')
    groups = lapply(seq_along(igraphL),  function(x) components(igraphL[[x]])$membership)
    
  }else(stop('please select from "rw",  "hcm", "km",  "pam",  "natrual" as method'))
  return(groups)
}

#' @title plot MCI barplots
#'
#' @description A barplot of MCI for all clusters in all states.
#' 
#' @param MCIl A list can be obtained through getMCI.
#' 
#' @param ylim A vector needed if the output barplots need to be on the same y scale.
#' 
#' @param nr The number of rows to plot.
#' 
#' @param nc The number of columns to plot,  default length(groups).
#' 
#' @param order A character vector of the order of the barplot.
#'  Default is NULL which uses the input list order.
#'  
#' @param minsize A non-negative numeric value of minimum size allowed for a cluster.
#' 
#' @param states A character of the state names to be shown on the plot. Default is NULL, 
#' assign this if you want to show all states including states without nodes.
#' 
#' @param title.size  Integer; the point size of the title. 
#'  This paramter is past to the graphical parameter \code{ps} to graphical function \code{title}. 
#'  What is meant by 'point size' is device-specific, but most devices mean a multiple of 1bp, that is 1/72 of an inch.  
#' 
#' 
#' @export
#' @return Return a barplot of MCI scores across states.
#' 
#' @references Chen L, Liu R, Liu Z, Li M & Aihara K (2012) 
#' Detecting early-warning signals for sudden deterioration of complex diseases by dynamical network biomarkers
#' Scientific Reports 2:342
#' 
#' @examples
#' test = list('state1' = matrix(sample(1:10, 6), 4, 3), 
#'   'state2' = matrix(sample(1:10, 6), 4, 3), 
#'   'state3' = matrix(sample(1:10, 6), 4, 3))
#' # assign colnames and rownames to the matrix
#' for(i in names(test)){
#'  colnames(test[[i]]) = 1:3
#'  row.names(test[[i]]) = c('g1', 'g2', 'g3', 'g4')
#' }
#'
#' cluster = list(c(1, 2, 2, 1),  c(1, 2, 3, 1),  c(2, 2, 1, 1))
#' names(cluster) = names(test)
#' for(i in names(cluster)){
#'  names(cluster[[i]]) = c('g1', 'g2', 'g3', 'g4')
#' }
#' membersL <- getMCI(cluster, test)
#' plotBar_MCI(membersL)
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

plotBar_MCI = function(MCIl, ylim = NULL, nr = 1, nc = NULL, 
                       order = NULL,  minsize = 3, states = NULL,
                       title.size = 30) ## update 08/28/2020 
{
  membersL = MCIl[[1]]
  MCI = MCIl[[2]]
  
  if(is.null(order)){
    loop = names(membersL)
  }else{
    if(any(!order %in% names(membersL)))
      stop('make sure all names provided in "order" is included in names of "countsL"')
    loop = order
  }
  if(!is.null(states)) loop = states
  if(is.null(nc)) nc = length(loop)
  par(mfrow = c(nr, nc))
  for(i in loop){
    if(! i %in% names(MCI)){
      mci = m = 0
    }else{
      mci = MCI[[i]]
      m = membersL[[i]]
      tmp = names(mci[is.na(mci)])
      if(length(tmp) != 0) m = m[!m %in% tmp]  # added in 10/27/2019
      nmembers = sapply(names(table(m)), function(x) length(m[m ==  x]))
      #cex = ifelse(length(m)>20, 0.7, 1)
      mci = mci[!is.na(mci)]
      
      if(!minsize<0 & minsize != 1){
        mci = mci[!nmembers<minsize]
        nmembers = nmembers[!nmembers<minsize]
      }else{
        warning('"minisize" need to be a non')
      }
      if(length(mci) ==  0) mci = 0
    }
    mci[is.na(mci)] = 0
    bar = barplot(mci, col = rainbow(length(mci),  alpha = 0.3), 
                  #legend = paste0('#', names(m), ' = ', sapply(m, nrow)), 
                  #main = paste0(i, ' (n = ',  max(m), ')'),  
                  main = '',  ## update 08/28/2020
                  ylab = 'MCI', 
                  xlab = 'modules', #args.legend = list(cex = cex)
                  ylim = ylim, 
                  cex.axis = 1.5,  
                  cex.names = 1.5, 
                  cex.main = 1.5,  
                  cex.lab = 1.5)
    title(main = paste0(i, '\n',  max(m), ' modules'), ps = title.size) ## update 08/28/2020
    if(all(mci != 0)) text(bar, mci, nmembers, cex = 1.5)
  }
}

#' @title Identifying the 'Biomodule'
#'
#' @description This function reports the 'biomodule', which is the module with
#'   the maximum Module Critical Index (MCI) scores for each state. Each state
#'   can have multiple modules (groups of subnetworks derived from the function
#'   \code{\link{getCluster_methods}}). This function runs over all states.
#'
#' @param membersL A list of vectors with unique sample ids as cluster names. The length
#'   of this list is equal to the number of states in the study. This can be the
#'   first element of the output from function \code{getMCI} or the output from
#'   \code{getCluster_methods}, see Examples for more detail.
#'   
#' @param MCIl A list of numeric vectors with unique cluster numbers as names.
#'   Each vector represents the MCI scores of that module. This can be the
#'   second element of the output from function \code{getMCI}.
#'   
#' @param minsize A numerical value of the minimum module size (the number of
#'   transcripts in a cluster) to output for downstream analysis.
#'
#' @return A nested list whose length is the length of the input object
#'   \code{membersL}.  Each internal list contains two objects: one object is
#'   the vector of biomodule IDs across states, and the other object is a list of
#'   transcript IDs (each defines the biomodule per state) across states.
#'   
#' @export
#' 
#' @examples
#' #1st option: get the input directly from getMCI function
#' test = list('state1' = matrix(sample(1:10, 6), 4, 3),  
#'      'state2' = matrix(sample(1:10, 6), 4, 3), 
#'      'state3' = matrix(sample(1:10, 6), 4, 3))
#'
#' # assign colnames and rownames to the matrix
#' for(i in names(test)){
#'   colnames(test[[i]]) = 1:3
#'     row.names(test[[i]]) = c('g1', 'g2', 'g3', 'g4')}
#'
#' cluster = list(c(1, 2, 2, 1), c(1, 2, 3, 1), c(2, 2, 1, 1))
#' names(cluster) = names(test)
#' for(i in names(cluster)){
#'   names(cluster[[i]]) = c('g1', 'g2', 'g3', 'g4')}
#'
#' membersL <- getMCI(cluster, test)
#' maxMCIms <- getMaxMCImember(membersL[[1]],  membersL[[2]],  min = 3)
#' #The same as
#' maxMCIms <- getMaxMCImember(cluster,  membersL[[2]],  min = 2)
#'
#'## case1: using 'rw' method by default
#'igraphL <- getNetwork(test,  fdr = 1)
#'cl <- getCluster_methods(igraphL)
#'## make sure every element in list cl is a \code{communities} object
#'sapply(cl, class)
#'##       state1        state2        state3
#'##"communities" "communities" "communities"
#'
#'## If there is(are) state(s) that is(are) empty which will not be a communities object(s),  
#'## please manually remove that state(s).
#'
#'cl = cl[which(sapply(cl, class) ==  'communities')]
#'
#'## and then run
#'library(igraph)
#'cluster = lapply(cl,  membership)
#'maxCIms <- getMaxMCImember(cluster,  membersL[[2]],  min = 2)
#'
#'## or run function 'getMCI' and use the 1st option
#'membersL <- getMCI(cl, test)
#'
#'## case2: using methods other than the default
#'cl <- getCluster_methods(test, method = "pam", cutoff = 2)
#'## check to make sure membersL[[2]] has values and run
#'maxCIms <- getMaxMCImember(cl,  membersL[[2]],  min = 2)
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}, Xinan Yang \email{xyang2@@uchicago.edu}

getMaxMCImember = function(membersL, MCIl, minsize = 1, n = 1)
{
  if(n<1 | class(n) != "numeric") stop('please provide a >= 1 numeric for n')
  if(is.null(names(membersL))) names(membersL) <- 1:length(membersL)
  n <- round(n)
  listn = names(membersL)
  
  if(minsize >= 1){
    minsize = minsize-1
    CIl = lapply(seq_along(membersL), function(x)
      ifelse(table(membersL[[x]])>minsize, MCIl[[x]], NA))
    module_keep = lapply(seq_along(membersL),  function(x)
      names(table(membersL[[x]])[table(membersL[[x]])>(minsize)]))  # corrected by xy from minsize-1 on 6/25/2020!
    membersL = lapply(seq_along(membersL), function(x)
      membersL[[x]][membersL[[x]] %in% module_keep[[x]]])
  } else {
    stop('please provide a minimum number of the cluster of interest,  
         which should be an integer that is larger than 0')
  }
  
  if(n>= 1) {
    idx = lapply(CIl, which.max)  # corrected on 8/28/2020
    maxCI = lapply(seq_along(idx), function(x) names(membersL[[x]][membersL[[x]] ==  idx[x]]))
    names(maxCI) = listn
    names(idx) = listn
    results <- list(idx = idx, members = maxCI)
  } 
  names(CIl) = listn
  if(n>1){
    for(j in 2:n){
      x <- unlist(lapply(idx, length))
      x <- names(x)[x>0]
      if(length(x)>0){
        for (i in x){
          CIl[[i]][idx[[i]]] <- NA  # mask the topest MCI score
          idx[[i]] = c(idx[[i]], which.max( CIl[[i]]))  # looking for the next maximum
        }
      }
      results[[j+1]] <- lapply(seq_along(idx), function(x) unlist(names(membersL[[x]][membersL[[x]] ==  idx[[x]][j]]))) 
      names(results[[j+1]]) = listn
      names(results)[j+1] <- paste0(j,'topest.members')
    }
  }
  
  results[['idx']] <- idx
  
  return(results)
}


#' @title  Get the cluster index and network nodes of biomodule
#'
#' @description This function retrieves the cluster index and network-node ids 
#' for the identified biomodule (that shows the maximum MCI score) at each state in the study.
#'
#' @param membersL A two-layer nested list of character or numeric values,  
#' any one out of the five elements output by the function \code{\link{getMCI}}.
#' 
#' @param idx A vector of integers that are cluster ids of the biomodule 
#' (the module with the highest MCI score) per state. 
#' This is the first element of the result from \code{\link{getMaxMCImember}}.
#' 
#' @return A list describing the biomodule of each state,  corresponding to one of the five elements 
#' (members,  MCI,  Sd,  PCC,  and PCCo) outputted by the function \code{\link{getMCI}}. 
#' The class of the vector depends on the class of the input parameter \code{membersL}.
#' 
#' @export
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @seealso \code{\link{getMCI}}
#' 
#' @examples
#' set.seed(2020)
#' test = list('state1' = matrix(rnorm(50,0,1), 10, 5), 
#'             'state2' = matrix(rnorm(30,0,3), 10, 3), 
#'             'state3' = matrix(rnorm(40,0,1), 10, 4))
#' 
#' ## Assign colnames and rownames to the matrix
#' for(i in 1:length(test)){
#'   colnames(test[[i]]) = paste0('Cell_',i,1:ncol(test[[i]]))
#'   row.names(test[[i]]) = paste0('g',1:10)
#' }
#' 
#' gene.cluster = list(rep(1:2, 5), c(rep(1:3,3),2), rep(1:5,2))
#' names(gene.cluster) = names(test)
#' for(i in names(gene.cluster)){
#'   names(gene.cluster[[i]]) = paste0('g',1:10)
#' }
#' 
#' membersL <- getMCI(gene.cluster, test)
#' names(membersL)
#' # [1] "members" "MCI"     "sd"      "PCC"     "PCCo" 
#' 
#' ## A list of index of interested gene.cluster IDs per state
#' idx = c(which.max(membersL[['MCI']][[1]]),
#'        which.max(membersL[['MCI']][[2]]),
#'        which.max(membersL[['MCI']][[3]]))

#' names(idx) = names(membersL[['sd']])
#' getMaxStats(membersL[['sd']], idx)
#' getMaxStats(membersL[['MCI']], idx)
#' 

getMaxStats = function(membersL, idx)
{
  if(any(is.null(names(idx)))| any(!names(idx) %in% names(membersL))) 
    stop('please make sure "idx" has names and all of its names is included in names of "membersL"')
  member_max = lapply(names(idx), function(x) membersL[[x]][idx[[x]]])
  names(member_max) = names(idx)
  member_max = member_max[lengths(member_max)>0]
  member_max = sapply(member_max, function(x) mean(x[[1]]))
  return(member_max)
}



#' @title Get the transcript ID and statistics for the n top MCI scores
#' 
#' @description This function generates a list of the n top MCI scores.
#' 
#' @param min A numerical value of the minimum number of transcripts in a cluster. A cutoff that determines
#' the smallest cluster in transcripts numbers in downstream analysis. 
#' 
#' @param n An integer determines how many modules are evaluated (the top n MCI scores would be evaluated
#' for downstream analysis). 
#' 
#' @param modulesL A list of integer named vectors. 
#' The length of this list is equal to the number of states in the study.
#' The names of a vector are the ids of modules (gene clusters) per state. 
#' This can be the first element of the output from the functions getMCI() or getCluster_methods(), 
#' see Examples for more detail.
#' 
#' @param membersL A two-layer nested list of characters or numbers, being any one out of the five outputting elements 
#' by the function getMCI(). This parameter is used to extract the stats (MCI, standard deviation, etc.), 
#' which is used to calculate the top MCI's. 
#' 
#' @param MCI1 A list of numeric vectors with unique cluster numbers as names. Each vector represents the 
#' MCI scores of that module. This can be the second element of the outputting of the function getMCI().
#' 
#' @return Returns a numeric vector of the n top MCI scores
#' 
#' @export
#' @author Yuxi Sun \email{ysun11@@uchicago.edu}
#' 
#' @examples
#' #' test = list('state1' = matrix(sample(1:10, 6), 4, 3),  
#'      'state2' = matrix(sample(1:10, 6), 4, 3), 
#'      'state3' = matrix(sample(1:10, 6), 4, 3))
#'
#' # assign colnames and rownames to the matrix
#' for(i in names(test)){
#'   colnames(test[[i]]) = 1:3
#'     row.names(test[[i]]) = c('g1', 'g2', 'g3', 'g4')}
#'
#' cluster = list(c(1, 2, 2, 1), c(1, 2, 3, 1), c(2, 2, 1, 1))
#' names(cluster) = names(test)
#' for(i in names(cluster)){
#'   names(cluster[[i]]) = c('g1', 'g2', 'g3', 'g4')}
#'
#' membersL <- getMCI(cluster, test)
#' topMCI <- getTopMCI(membersL[[1]],  membersL[[2]],  membersL[[2]],
#' min = 3, n = 1)

getTopMCI = function(modulesL, MCI1, membersL, min, n = 1)
{
  maxMCIms <- getMaxMCImember(modulesL, MCI1, min) 
  topMCI = getMaxStats(membersL,maxMCIms[[1]])
  topMCI = topMCI[order(topMCI, decreasing = TRUE)]
  topMCI = topMCI[1:n] 
  return(topMCI)
}

#' @title Plot the Maximized MCI per State
#'
#' @description This function generates a line plot over multiple states with the maximum MCI score per state. 
#' The module size (i.e., number of network nodes) is specified at each state in parentheses.
#'
#' @param maxMCIms A list of 2 elements. The 1st element is an integer vector of module ids whose names are the state names. 
#' The 2nd element is a list of character vectors per state. The vectors are network nodes (e.g. transcript ids). 
#' This parameter can be obtained by running function \code{\link{getMaxMCImember}}.
#' 
#' @param MCIl A list of numeric vectors whose names are unique cluster ids. 
#' Each vector represents the MCI scores of modules in a state. 
#' This can be the second element of the output from the function \code{\link{getMCI}}.
#' 
#' @param las Numeric in {0, 1, 2, 3}; the style of axis labels. Default is 0, meaning labels are parallel. 
#' See \code{\link{getMCI}} for more detail.
#' 
#' @param order A vector of state names in the customized order to be plotted, setting to NULL by default.
#' 
#' @param states A character vector of state names that will be shown on the plot, setting to NULL by default.
#' Assign this if you want to show all states, including states with no resulted modules. 
#' This parameter will overwrite the parameter 'order'.
#' 
#' @return Returns a line plot of maximum MCI scores across the states
#' 
#' @export
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' 
#' @examples
#' maxMCIms = list(c(state1 = 1,  state2 = 2,  state3 = 1), 
#'    c(list(state1 = c('g1', 'g2', 'g3'), 
#'    state2 = c('g3', 'g5'), 
#'    state3 = c('g2', 'g6'))))
#'    
#' MCIl = list(state1 = c('1' = 8.84, '2' = 6.4),  
#'    state2 = c('1' = NA, '2' = 9.5, '3' = NA),    
#'    state3 = c('1' = 2.3, '2' = 1.4))
#'    
#' plotMaxMCI(maxMCIms, MCIl)
#'

plotMaxMCI = function(maxMCIms,  MCIl,  las = 0,  order = NULL,  states = NULL)
{
  if(any(is.null(names(maxMCIms[[1]]))) | any(is.null(names(maxMCIms[[2]])))) 
    stop('Please give names for the 1st and 2nd element of the "maxMCIms" as well as "MCIl"')
  if(is.null(order)){
    CI = sapply(names(maxMCIms[[1]]),  function(x) MCIl[[x]][maxMCIms[[1]][[x]]])
    ln = names(maxMCIms[[1]])
    names(CI) = ln
  }else{
    if(any(!order %in% names(maxMCIms[[2]])))
      stop('make sure all names in "order" are in names of the 2nd element of "maxMCIms"')
    if(any(!names(maxMCIms[[2]]) %in% order))
      warning('not every state in "simulation" is plotted,  make sure "order" is complete')
    CI = sapply(order,  function(x) MCIl[[x]][maxMCIms[[1]][[x]]])
    ln = order
  }
  if(any(is(CI) ==  'list')){
    warning('changing NA CI score(s) to 0')
    idx = sapply(CI, function(x) length(x) ==  0)
    CI[idx] = 0
    CI = do.call(c, CI)
    names(CI) = ln
  }
  if(!is.null(states)){
    CI = CI[states]
    CI[is.na(CI)] = 0
    ln = names(CI) = states
  }
  matplot(CI, type = 'l', ylab = 'MCI(m|r)', axes = FALSE)
  len = sapply(ln, function(x) length(maxMCIms[[2]][[x]]))
  len[is.na(len)] = 0
  names(len) = ln
  text(seq_along(CI), CI+0.01, paste0('(', len, ')'))
  
  axis(2)
  axis(side = 1, at = seq_along(CI), labels = ln, las = las)
}

#' Obtain the identified BioTiP and its length
#' @description Obtain the identified CTS transcripts.
#' 
#' @param maxMCI A list of numeric vectors,  whose lengths are the numbers of system's states. 
#' This gives the maximum MCI score of each state,  and it can be obtained from the output of \code{\link{getMaxStats}}. 
#' Names of the list need to be included in names of \code{maxMCIms}.
#' 
#' @param maxMCIms A list of character vectors whose length is the system's states studied. 
#' The vectors records constructed gene-gene coexpression network nodes (e.g. transcript ids). 
#' This parameter is the second element of the output of the function \code{\link{getMaxMCImember}}.
#' 
#' @return A list of character vectors, in which the elements are the unique IDs of the network nodes of the BioTiP.
#' 
#' @export
#' 
#' @examples
#' maxMCI <- c(a = 2.56,  b = 8.52,  c = 2.36,  d = 4.81,  e = 5.26)
#' maxMCIms <- list(a = c("A100",  "A293",  "C403"),  
#'                  b = c("B853",  "D826",  "A406"),  
#'                  c = c("J198",  "D103",  "B105"),  
#'                  d = c("K529",  "D385",  "E358"),  
#'                  e = c("J019",  "U926",  "N824"))
#' identical(names(maxMCI),  names(maxMCIms))
#' # TRUE
#' getCTS(maxMCI,  maxMCIms)
#' # "Length: 3"
#' # "B853" "D826" "A406"
#' 
#' @author Antonio Feliciano y Pleyto, Zhezhen Wang \email{zhezhen@@uchicago.edu},
#' Yuxi Sun \email{ysun11@@uchicago.edu}

getCTS <- function(maxMCI,  maxMCIms) 
{
  if (is.null(names(maxMCI))) {
    stop("No names for maxMCI. Please provide names.")
  }
  if (is.null(names(maxMCIms))) {
    stop("No names for maxMCIms. Please provide names.")
  }
  if (!all(names(maxMCI) %in% names(maxMCIms))) {
    stop("Names of maxMCI has to be in maxMCIms.")
  }
  CTS_list = vector(mode = "list", length = length(maxMCI))
  for (i in 1:length(maxMCI)) {
    CTS_list[[i]] <- maxMCIms[[names(maxMCI)[i]]]
    message(paste0("Length: ",  length(CTS_list[[i]])))
  }
  names(CTS_list) <- names(maxMCI)   ## added by Holly 06232020
  return(CTS_list)  
}

#' @title plot a line plot of Ic scores for each state.
#'
#' @description plot a line plot with Ic score for each state
#'
#' @param Ic A vector with names of states. If order is not assigned,  
#' then plot by the order of this vector.
#' 
#' @param las Numeric in {0, 1, 2, 3}; the style of axis labels. 
#' Default is 0,  meaning labels are parallel. 
#' (link to http://127.0.0.1:21580/library/graphics/html/par.html).
#' 
#' @param order A vector of state names in the customized order to be plotted, setting to NULL by default.
#' 
#' @param ylab titles y axes, as in plot.
#' 
#' @param col vector of colors. Colors are used cyclically.
#' 
#' @param main A character vector. The title of the plot. Default is NULL.
#' 
#' @param add logical. If TRUE,  plots are added to the current one. 
#' This is inherited from \link[graphics]{matplot}.
#' 
#' @param ylim An integer vector of length 2. Default is NULL.
#' 
#' @param lty An vector of line types. This is also inherited from \link[graphics]{matplot}.
#' 
#' @param lwd Anineger of line widths. This is also inherited from \link[graphics]{matplot}.
#' 
#' @export
#' 
#' @return Return a line plot of Ic score across states.
#' 
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' 
#' @examples
#' Ic = c('state3' = 3.4,  'state1' = 5.6,  'state2' = 2)
#' plotIc(Ic, order = c('state1', 'state2', 'state3'))

plotIc = function(Ic,  las = 0,  order = NULL,  ylab = "Ic",  col = "black",  
                  main = NULL,  add = FALSE,  ylim = NULL, lty = 1:5,  lwd = 1)
{
  if(!is.null(order)){
    if(any(!order %in% names(Ic))) stop('make sure "Ic" is named using names in "order"')
    if(any(!names(Ic) %in% order)) warning('not every state in "Ic" is plotted,  make sure "order" is complete')
    Ic = Ic[order]
  }
  matplot(Ic,  type = "l",  ylab = ylab,  axes = FALSE,  col = col,  main = main,  add = add,  ylim = ylim, lty = lty, lwd = lwd)
  axis(2)
  stages = names(Ic)
  axis(side = 1, at = seq_along(Ic), labels = stages, las = las)
}

#######################################################################
################# updated functions,  02/24/2020 #######################
#######################################################################

#' @title Building Networks of Nodes
#'
#' @description This function builds one correlation network for each state
#'   (sample group) and runs across all states. The network nodes are defined by
#'   the context of the input dataset. For transcriptomic network analysis, 
#'   network nodes can be the expressed transcript IDs and network links can be
#'   the correlation coefficients. Using the Pearson Correlation Coefficient
#'   (PCC) analysis, this function assembles a correlation network of nodes
#'   (e.g., co-expressed transcripts) for each state using the R package igraph.

#' @param optimal A list of x numeric data frames, where x is the number of
#'   states studied. Each data frame consists of loci with high standard
#'   deviations. This object can be obtained through \code{sd_selection}
#'   function.
#'   
#' @param fdr A numeric cutoff value for a Pearson Correlation Coefficient
#'   (PCC) analysis. Default is 0.05. Transcripts are linked into a network if
#'   their correlations meet this PCC-significance criterion.
#'
#' @return A list of igraph objects whose length is the length of the input
#'   object \code{optimal}. Each object is a network of correlated nodes whose
#'   PCCs meet the significant criteria based on the false discovery rate (FDR)
#'   control. The length of the list is the number of states with PCC networks.
#'   If no PCC meets the significant criteria in a state, the state
#'   will be deleted from the output.
#'   
#' @export
#' @import stringr psych
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}
#'
#' @examples
#' test = list('state1' = matrix(sample(1:10, 6), 2, 3), 
#'  'state2' = matrix(sample(1:10, 6), 2, 3), 
#'  'state3' = matrix(sample(1:10, 6), 2, 3))
#'
#' for(i in names(test)){
#'   colnames(test[[i]]) = 1:3
#'   row.names(test[[i]]) = 1:2}
#'
#' igraphL <- getNetwork(test,  fdr = 1)
#' #[1] "state1:2 nodes"
#' #[1] "state2:2 nodes"
#' #[1] "state3:2 nodes

getNetwork = function(optimal, fdr = 0.05)
{
  #  require(stringr)
  #  require(psych)
  #  require(igraph)
  
  rL = lapply(optimal, function(x) corr.test(t(x), adjust = 'fdr', ci = FALSE)$r)
  names(rL) = names(optimal)
  pL = lapply(optimal, function(x) corr.test(t(x), adjust = 'fdr', ci = FALSE)$p)
  if(is.null(names(rL))) stop('give names to the input list')
  
  igraphL = list()
  for(i in names(rL)){
    test = rL[[i]]
    test.p = pL[[i]]
    # temprally replace '.' in a gene name with '+" ## note # 1
    row.names(test) = gsub('.', '+', row.names(test), perl = FALSE, fixed = TRUE)   # modified on 2/24/2020 
    row.names(test.p ) = gsub('.', '+', row.names(test.p), perl = FALSE, fixed = TRUE)  # modified on 2/24/2020 
    test[lower.tri(test, diag = TRUE)] = NA
    #test.p[lower.tri(test, diag = TRUE)] = 1
    tmp = lapply(1:nrow(test), function(x) test[x, test.p[x, ]<fdr])
    tmp_name = lapply(1:nrow(test), function(x) which(test.p[x, ]<fdr))
    idx = which(lengths(tmp_name) ==  1)
    for(j in idx){
      names(tmp[[j]]) = names(tmp_name[[j]])
    }
    names(tmp) = row.names(test)
    edges = stack(do.call(c, tmp))
    edges = subset(edges,  !is.na(edges$values))
    #tmp2 = subset(edges, grepl('[.][1-9, A-z][.]', ind))
    #if(nrow(tmp2)!= 0){
    #  tmp2$node1 = paste0(str_split_fixed(tmp2$ind, '\\.', 3)[, 1], '.', str_split_fixed(tmp2$ind, '\\.', 3)[, 2])
    #  tmp2$node2 = str_split_fixed(tmp2$ind, '\\.', 3)[, 3]
    #}
    #edges = subset(edges, !grepl('\\.[1-9, A-z]\\.', ind))
    edges$node1 = str_split_fixed(edges$ind, '\\.', 2)[, 1]
    edges$node2 = str_split_fixed(edges$ind, '\\.', 2)[, 2]
    #edges = rbind(edges, tmp2)
    # return back the '.' in a gene name from temprally used '+" ## note # 1
    edges$node1 = gsub('+', '.', edges$node1, fixed = TRUE)  # modified on 2/24/2020
    edges$node2 = gsub('+', '.', edges$node2, fixed = TRUE)  # modified on 2/24/2020
    edges = edges[, c('node1', 'node2', 'values')]
    edges$weight = abs(edges$values) # added in 1/8/2019
    #colnames(edges) = c('node1', 'node2', 'weight') # added in 12/18/2018
    
    nodes = data.frame(unique(c(edges$node1, edges$node2)))
    message(paste0(i, ':', nrow(nodes), ' nodes')) #[1] 48    1
    routes_igraph <- graph_from_data_frame(d = edges,  vertices = nodes,  directed = FALSE)
    igraphL[[i]] = routes_igraph
  }
  return(igraphL)
}


#######################################################################
################# updated functions,  02/14/2020 #######################
#######################################################################


#' @title Calculating (and plot) random Ic scores (Mojtahedi et al. 2016) based on shuffling sample labelling.
#'
#' @description Run \code{B} times of sample-label shuffling to calculate the Ic score,  
#' where x should be the same as the length of identified BioTiP and B is self-defined.
#'
#' @param counts A numeric matrix or data frame. The rows and columns 
#' represent unique transcript IDs (geneID) and sample names,  respectively.
#' 
#' @param sampleNo An integer of sample size at the tipping-point state.
#' 
#' @param Ic A numeric value. Ic score of identified CTS (gene-set), useful when \code{plot} is TRUE.
#' Default is NULL.
#' 
#' @param genes A character vector of identified CTS gene unique ids.
#' 
#' @param B An integer indicating the number of times to run this simulation, default 1000.
#' 
#' @inheritParams getIc 
#' 
#' @param ylim An integer vector of length 2. Default is NULL.
#' 
#' @param main A character vector. The title of the plot. Default is NULL.
#' 
#' @param plot A Bollen value indicating whether a density plot will be plotted.
#' 
#' @export
#' 
#' @return A vector of \code{B} values of BioTIP (or Ic) scores calculated for the state of interest.
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}
#' 
#' @examples
#' counts = matrix(sample(1:100, 27), 3, 9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1', 'loci2', 'loci3')
#' CTS = c('loci1', 'loci2')
#' randomS <- simulation_Ic_sample(counts,  sampleNo = 3,  Ic = 3.4,  genes = CTS,  B = 3, 
#'                        fun = 'BioTIP', plot = TRUE)
#' dim(randomS)

## updated from the old function plot_simulations_sample()    2/19/2020
simulation_Ic_sample = function(counts,  sampleNo,  Ic = NULL,  genes, B = 1000, 
                                ylim = NULL, main = 'simulation of samples',
                                fun = c("cor", "BioTIP"),   
                                shrink = TRUE,  
                                use = c("everything",  "all.obs",  "complete.obs",  "na.or.complete",  "pairwise.complete.obs"),
                                output = c('Ic', 'PCCg', 'PCCs'),
                                plot = FALSE,
                                PCC_sample.target = 1)  ## 12/02/2020
{
  output <- match.arg(output)
  fun <- match.arg(fun)
  use <- match.arg(use)
  #PCC_sample.target = 'average'
  # begin ## 12/02/2020
  # PCC_sample.target = match.arg(PCC_sample.target)
  if (class(PCC_sample.target) ==  'numeric') if((PCC_sample.target < 0) | (PCC_sample.target > 1)) 
  { 
    stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero', 'average', 'half'")
  } else if(class(PCC_sample.target) ==  'character') if(!PCC_sample.target %in% c('none', 'zero',  'average', 'half')) {
    stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero', 'average', 'half'")
  }  
  # end ## 12/02/2020
  PCC_gene.target = 'zero'
  
  # random select sampleNo cells regardless its labelling for state
  sampleL = lapply(1:B,  function(x) sample(colnames(counts), sampleNo))
  tmp = sapply(1:B,  function(x) 
    getIc(counts, sampleL = sampleL[x], genes = genes, output = output,
          fun = fun,  
          shrink = shrink,
          use = use,
          PCC_sample.target = PCC_sample.target))
  
  if(plot) {
    p_v = length(tmp[tmp>Ic])/B
    den = density(tmp)
    xmin = min(Ic, den$x)
    xmax = max(Ic, den$x)
    plot(den, main = main, xlim = c(xmin, xmax), ylim = ylim)
    abline(v = Ic, col = 'red', lty = 2)
    x = max(den$x) - 0.2*diff(range(den$x))
    if(p_v ==  0) p_v = paste('<', 1/B)
    text(x, max(den$y)-0.05, paste('P = ', p_v))
  }
  return(tmp)
}


#' @title Get Index for Critical transition (Ic score)
#'
#' @description Retrieve Ic scores (Pearson correlation of genes / Pearson correlation of samples) 
#' for the identified critical transition state
#' 
#' @param counts A  numeric matrix or data frame. The rows and columns represent unique transcript IDs (geneID)
#'  and sample names,  respectively.
#' 
#' @param sampleL A list of vectors, whose length is the number of states. Each vector gives the sample names in a state. 
#' Note that the vector s (sample names) has to be among the column names of the R object 'df'.
#' 
#' @param genes A character vector consisting of unique CTS gene ids. This can be obtained from \code{\link{getMaxMCImember}}
#' 
#' @param output A string. Please select from 'Ic',  'PCCg',  or 'PCCs'. Uses 'Ic' by default.
#' 'PCCg' is the PCC between genes (numerator) and 'PCCs' is PCC between samples (denominator)
#' 
#' @param fun An optional character string indicating the R functon to calculate correlations 
#' for all possible pairs of columns of a matrix. 
#' When using "BioTIP",  The method is modified to ignore missing values, analogous to how
#' \code{cor(X,  use = "pairwise.complete.obs")} works.  
#' Note that the "BioTIP" option only function together with shrink = TRUE.
#' 
#' @param shrink A flag specifying whether to shrink the matrix of gene-gene correlation or not. 
#' This appraoch uses the method outlined by Schafer and Strimmer in 
#' "A Shrinkage Approach to Large-Scale Covariance Matrix Estimation 
#' and Implications for Functional Genomics" (2005). Here, we shrink between-gene correlations 
#' towards 0 due to the low global gene expressional dependence in a stable state 
#' Comparing to fun = 'cor', the 'BioTIP' method without shinkage is modified 
#' to ignore missing values, analogous to how \code{cor(X, use = "pairwise.complete.obs")} works. 
#' For between-sample correlation matrix, we shrink 
#' towards the average correlation to reflect the similar gene-expression profiles in a stable state. 
#' 
#' @param use An optional character string,  when fun ==  "cor",  it gives a method 
#' for computing covariances in the presence of missing values. 
#' This must be (an abbreviation of) one of the strings "everything", "all.obs", 
#' "complete.obs", "na.or.complete", or "pairwise.complete.obs". 
#' 
#' @param PCC_sample.target A numeric choose between [0,1] or 
#' a character choose among ('none, 'average',  'zero', 'half'),  
#' indicating whether to turn off shrinkage, or to shrink PCC towards their empirical common average,  
#' zero or 0.5 (for sample-sample correlations).
#' 
#' @return A list of numeric values,  whose length and names are inherited from \code{sampleL}
#' @export
#' 
#' @references Schafer and Strimmer (2005) "A Shrinkage Approach to Large-Scale 
#' Covariance Matrix Estimation and Implications for Functional Genomics"
#' 
#' @references M. Mojtahedi et al., Cell Fate Decision as High-Dimensional Critical State Transition. 
#' PLoS Biol 14,  e2000640 (2016).
#' 
#' @examples
#' counts = matrix(sample(1:100, 27),  3,  9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1', 'loci2', 'loci3')
#' cli = cbind(1:9, rep(c('state1', 'state2', 'state3'), each = 3))
#' colnames(cli) = c('samples', 'group')
#' samplesL <- split(cli[, 1],  f = cli[, 'group'])
#' CTS = c('loci1', 'loci2')
#' 
#' ## Comparing the results with an estiamted correlation matrix with that without estimation.
#' Ic = getIc(counts,  samplesL,  CTS, fun = 'cor')
#' Ic.2 = getIc(counts,  samplesL,  CTS, fun = 'BioTIP', shrink = FALSE)
#' BioTIP = getIc(counts,  samplesL,  CTS, fun = 'BioTIP')
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}

getIc <- function(counts,  sampleL,  genes,  output = c('Ic', 'PCCg', 'PCCs'), 
                  fun = c("cor", "BioTIP"), 
                  shrink = TRUE, 
                  use = c("everything",  "all.obs",  "complete.obs",  "na.or.complete",  "pairwise.complete.obs"),
                  PCC_sample.target = 1)    ## 12/02/2020
{
  if (class(genes)!='character') stop("genes have to be a character of gene symbols, i.e. 
                                      genes have to be a subset of row.names(counts)")
  output <- match.arg(output)
  fun <- match.arg(fun)
  use <- match.arg(use)
  #PCC_sample.target = 'average'
  # begin ## 12/02/2020
  # PCC_sample.target = match.arg(PCC_sample.target)
  if (class(PCC_sample.target) ==  'numeric') if((PCC_sample.target < 0) | (PCC_sample.target > 1)) 
  { 
    stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
  } else if(class(PCC_sample.target) ==  'character') if(!PCC_sample.target %in% c('none', 'zero',  'average', 'half')) {
    stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
  }  
  # end ## 12/02/2020
  
  PCC_gene.target = 'zero'
  
  subsetC = subset(counts, row.names(counts) %in% genes)
  subsetC = lapply(sampleL,function(x) subsetC[,as.character(x)])
  # if(fun ==  "BioTIP" & PCC_sample.target ==  'none') ## 12/02/2020 Now, this is a case we want 
  #   warning('You are not really calling BioTIP function without a proper setting of PCC_sample.target !') ## 12/02/2020
  if(fun ==  "BioTIP" & PCC_gene.target ==  'none') 
    warning('You are not really calling BioTIP function without a proper setting of PCC_gene.target !')
  
  # for "pairwise.complete.obs", remove those results of two or less pairs
  
  if (fun ==  "BioTIP") {
    PCCg = lapply(subsetC, function(x) avg.cor.shrink(x,
                                                      MARGIN = 1,
                                                      shrink = shrink,
                                                      abs = TRUE,
                                                      target = PCC_gene.target))
    PCCg = unlist(PCCg) 
  } else {
    PCCg = lapply(subsetC, function(x) abs(cor(t(x), use = use))) 
    for (i in seq_along(PCCg)) PCCg[[i]][upper.tri(PCCg[[i]], 
                                                   diag = FALSE)]   ## updated 02/17/20
    PCCg = sapply(PCCg, function(x) mean(x, na.rm = TRUE))
  }
  #if (fun ==  "BioTIP" & PCC_sample.target!= 'none') {# 1/4/2021
  if (fun ==  "BioTIP" ) { # 1/4/2021, allowing NA values to be included in subsetC when fun='BioTIP'
    PCCs = lapply(subsetC, function(x) avg.cor.shrink(x, 
                                                      MARGIN = 2, 
                                                      shrink = shrink,
                                                      abs = FALSE, 
                                                      target = PCC_sample.target ))   
    PCCs = unlist(PCCs)                      
  } else {
    PCCs = lapply(subsetC, function(x) cor(x, use = use))
    for (i in seq_along(PCCs)) PCCs[[i]][upper.tri(PCCs[[i]], 
                                                   diag = FALSE)]    ## updated 02/17/20
    PCCs = sapply(PCCs, function(x) mean(x, na.rm = TRUE))
  }
  toplot = PCCg/PCCs
  names(toplot) = names(PCCg) = names(PCCs) = names(sampleL)
  if (output ==  "Ic") {
    return(toplot)
  }
  else if (output ==  "PCCg") {
    return(PCCg)
  }
  else if (output ==  "PCCs") {
    return(PCCs)
  }
}



#' @title Calculating random Index of Critical transition (Ic scores) for randomly-selected genes
#'
#' @description Simulating Ic scores for \code{x} randomly selected samples, where x should be the same 
#' as the length of the identified critical-transition signal (CTS) (e.g.,  number of genes) and \code{B} is self-defined running times.
#'
#' @param obs.x An integer, length of identified CTS.
#' 
#' @param sampleL A list of vectors,  whose length is the number of states. 
#' Each vector gives the sample names in a state. Note that the vector s (sample names) 
#' has to be among the column names of the R object 'df'.
#' 
#' @param counts A numeric matrix or dataframe in which columns are samples and rows are transcripts.
#' Each row needs to have a unique row name (i.e. transcript ID).
#' 
#' @param B An integer, setting the permutation with \code{B} runs. Default is 1000.
#' 
#' @inheritParams getIc
#'   
#' @param use An optional character string,  when fun ==  "cor", it gives a method 
#' for computing covariances in the presence of missing values. 
#' This must be (an abbreviation of) one of the strings "everything",  "all.obs",  
#' "complete.obs",  "na.or.complete",  or "pairwise.complete.obs". 
#' 
#' @return A matrix of \code{y} rows and \code{B} columns where \code{y} 
#' is the length of \code{sampleL} and \code{B} is self-defined. Each column is a set of Ic scores calculated for each state
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}
#' @export
#' 
#' @examples
#' counts = matrix(sample(1:100, 27), 3, 9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1', 'loci2', 'loci3')
#' cli = cbind(1:9, rep(c('state1', 'state2', 'state3'), each = 3))
#' colnames(cli) = c('samples', 'group')
#' samplesL <- split(cli[, 1], f = cli[, 'group'])
#' simulation_Ic(2, samplesL, counts, B = 3, fun = "BioTIP", shrink = TRUE)

simulation_Ic <- function(obs.x,  sampleL,  counts,  B = 1000,  fun = c("cor", "BioTIP"),   
                          shrink = TRUE, 
                          use = c("everything",  "all.obs",  "complete.obs",  "na.or.complete",  "pairwise.complete.obs"),
                          output = c('Ic', 'PCCg', 'PCCs'),
                          PCC_sample.target = 1) ## 12/02/2020
{
  fun <- match.arg(fun)
  use <- match.arg(use)
  
  #PCC_sample.target = 'average'
  # begin ## 12/02/2020
  # PCC_sample.target = match.arg(PCC_sample.target)
  if (class(PCC_sample.target) ==  'numeric') if((PCC_sample.target < 0) | (PCC_sample.target > 1)) 
  { 
    stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
  } else if(class(PCC_sample.target) ==  'character') if(!PCC_sample.target %in% c('none', 'zero',  'average', 'half')) {
    stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
  }  
  # end ## 12/02/2020
  
  PCC_gene.target = 'zero'
  output <- match.arg(output)
  
  #  set.seed(2020)
  random = sapply(1:B,  function(x) sample(row.names(counts), obs.x))
  
  # create progress bar
  pb <- txtProgressBar(min = 0,  max = B,  style = 3)
  m <- matrix(nrow = length(sampleL),  ncol = B)
  for(i in 1:B)
  {
    setTxtProgressBar(pb,  i)
    m[, i] <- getIc(counts, sampleL = sampleL, genes = random[, i], output = output,  ## updated 02/17/20
                    fun = fun,  
                    shrink = shrink, 
                    use = use,
                    PCC_sample.target = PCC_sample.target)
    Sys.sleep(0.01)
    if(i ==  B) cat("Done!\n")
  }
  # m = sapply(1:B,  function(x) getIc(counts, sampleL, random[, x], output = 'Ic',  
  #                                    fun = fun,  shrink = TRUE, 
  #                                    PCC_sample.target = PCC_sample.target,  
  #                                    PCC_gene.target = PCC_gene.target, 
  #                                    use = use))
  
  row.names(m) = names(sampleL)
  return(m)
}



#' @title Line or boxplot of an observed and its simulated scores
#'
#' @description Generate a line (or box) plot of Ic score and simulated Ic scores, 
#' with three horizontal lines: the min,  max and 2*(max-min) value of the state of interests, or all values.
#'
#' @inheritParams plotIc 
#' 
#' @param simulation A numeric matrix of Ic scores in which rows are states and columns are simulation runs. 
#' It can be obtained from \code{\link{simulation_Ic}}
#' 
#' @param order Characters of names of Ic to be plotted in a desired \code{order}. Default is NULL.
#' 
#' @param fun A character choose between ('matplot', 'boxplot'), indicating plot type.
#' 
#' @param which2point A character (or integer) which state's values were used to set up the three horizontal lines. 
#' by default is NULL,  indicating the values of all states will be used.
#' 
#' @export
#' 
#' @return Return a plot of the observed Ic (red) and simulated Ic (grey) scores per state.
#' 
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}
#' 
#' @examples
#' sim = matrix(sample(1:10, 9), 3, 3)
#' row.names(sim) = paste0('state', 1:3)
#' Ic = c('state1' = 3.4, 'state2' = 5.6, 'state3' = 2)
#' plot_Ic_Simulation(Ic, sim)

plot_Ic_Simulation <- function (Ic,  simulation,  las = 0,  ylim = NULL,  
                                order = NULL,  main = NULL,  
                                ylab = "Ic",  fun = c('matplot', 'boxplot'), 
                                which2point = NULL) 
{
  fun <- match.arg(fun)
  if (any(is.null(names(Ic)))) 
    stop("Please provide name for vector \"Ic\" ")
  if (any(is.null(rownames(simulation)))) 
    stop("Please provide rowname for vectors of \"simulation\" ")
  if(length(Ic)!= nrow(simulation)) 
    stop("Please provide the same length of \"Ic\" and vectors of \"simulation\" ")
  if (!identical(names(Ic),  row.names(simulation))) 
    Ic = Ic[match(row.names(simulation),  names(Ic))]
  if(!is.null(which2point)) {
    if (! (which2point %in% rownames(simulation) | which2point %in% 1:nrow(simulation))) 
      stop("which2point must be a state name of integer indicating the state of interested.")
  }
  
  if (fun ==  'matplot') {
    toplot = cbind(simulation,  Ic)
  } else { toplot = simulation }
  if (!is.null(order)) {
    if (any(!names(Ic) %in% order)) 
      warning("not all states in Ic is plotted")
    if (any(!order %in% names(Ic))) 
      stop("make sure \"Ic\" is named using names in \"order\"")
    toplot = toplot[order,  ]
  }
  if(is.null(ylim)) {
    if(is.null(which2point)) {
      ylim = c(min(c(Ic, simulation)),  max(Ic,  2*(max(simulation)-min(simulation))))
    } else {
      ylim = c(min(c(Ic, simulation)),  
             max(max(Ic[which2point],  
                     2*(max(simulation[which2point, ])-min(simulation[which2point, ])),
                     simulation)))
    }
  }
  if (fun ==  'matplot') {
    matplot(toplot,  type = "l",  
            col = c(rep("grey",  ncol(toplot) - 1),  "red"),  
            lty = 1,  ylab = ylab,  
            axes = FALSE,  ylim = ylim,  main = main)
  } else {
    boxplot(t(toplot),  col = c(rep("grey",  ncol(toplot) - 1),  "red"),  
            ylab = ylab,  
            axes = FALSE,  ylim = ylim,  main = main)
    points(Ic,  col = "red",  type = 'b')
    x <- lapply(Ic,  function(x) table(toplot>x))
    y <- unlist(lapply(x,  function(X) X[2]/sum(X)))
    if(any(is.na(y))) { y[which(is.na(y))] = 0 }
    sig <- which(y<0.05)
    if(length(sig)>0) mtext( round(y[sig], 3),  
                             line = -5,  at = (1:length(Ic))[sig])
  }          
  
  axis(2)
  stages = row.names(toplot)
  axis(side = 1,  at = seq_along(stages),  labels = stages,  las = las)
  
  if(is.null(which2point)) {
    abline(h = min(simulation),  col = "grey",  lty = 3)
    abline(h = max(simulation),  col = "grey",  lty = 3)
    abline(h = min(simulation)+ 2*(max(simulation)-min(simulation)),  
           col = "grey",  lty = 2)
  }  else {
    abline(h = min(simulation[which2point, ]),  col = "grey",  lty = 3)
    abline(h = max(simulation[which2point, ]),  col = "grey",  lty = 3)
    abline(h = min(simulation[which2point, ])+ 
             2*(max(simulation[which2point, ])-min(simulation[which2point, ])),  
           col = "grey",  lty = 2)
  }
  
  
}




#' @title Calculating MCI Scores
#'
#' @description This function calculates a module critical index (MCI) score for
#'   each module per state within a dataset. Each module is a cluster of
#'   transcripts generated from the function \code{\link{getCluster_methods}}.
#'   Note that a dataset should contain three or more states (samples in
#'   groups).
#'
#' @param groups A list of elements whose length is the member of states. The
#'   elements could be either be vectors or \code{communities} object of the
#'   R package \code{\link{igraph}}. If a vector, it is the output of the function
#'   \code{getCluster_methods}. The names of each vector are the pre-selected
#'   transcript IDs generated by the function \code{\link{sd_selection}}. Each
#'   vector, whose length is the number of pre-selected transcripts in a state, 
#'   contains the module IDs. If a \code{communities} object, it can be obtained
#'   by \code{getCluster_methods} using the "rw" method. It is also an output of
#'   the function \code{\link{sd_selection}}.
#'   
#' @param countsL A list of x numeric count matrices or x data frame, where x is
#'   the number of states.
#'   
#' @param adjust.size A Boolean value indicating if MCI score should be adjusted
#'   by module size (the number of transcripts in the module) or not. Default FALSE.
#'   This parameter is not recommended for fun = BioTIP.
#'   
#' @param fun A character chosen between ("cor", "BioTIP"), indicating where an adjusted 
#'   correlation matrix will be used to calculate the MCI score.   
#'
#' @param df NULL or a numeric matrix or data frame, where rows and columns represent
#'   unique transcript IDs (geneID) and sample names, respectively. 
#'   Used only when \code{fun = 'BioTIP'}. 
#'   By default is NULL, estimating the correlation among selected genes. 
#'   Otherwise, estimating the correlation among all genes in the df, ensuring cross-state comparison.    
#'
#' @return A list of five elements (members,  MCI,  Sd,  PCC,  and PCCo). Each of
#'   element is a two-layer nested list whose length is the length of the input
#'   object \code{groups}. Each internal nested list is structured according to
#'   the number of modules identified in that state.
#' * members: vectors of unique ids
#' * MCI: the MCI score
#' * sd: standard deviation
#' * PCC: Mean of pairwise Pearson Correlation Coefficient calculated among the
#' loci in a module.
#' * PCCo: Mean of pairwise Pearson Correlation Coefficient calculated between
#' the loci in a module and the loci outside that module but inside the same
#' state.
#' @export
#'
#' @examples
#' test = list('state1' = matrix(sample(1:10, 6), 4, 3), 'state2' = 
#' matrix(sample(1:10, 6), 4, 3), 'state3' = matrix(sample(1:10, 6), 4, 3))
#'
#' ## Assign colnames and rownames to the matrix
#' for(i in names(test)){
#'    colnames(test[[i]]) = 1:3
#'    row.names(test[[i]]) = c('g1', 'g2', 'g3', 'g4')}
#'
#' cluster = list(c(1, 2, 2, 1), c(1, 2, 3, 1), c(2, 2, 1, 1))
#' names(cluster) = names(test)
#' for(i in names(cluster)){
#'    names(cluster[[i]]) = c('g1', 'g2', 'g3', 'g4')}
#'
#' membersL <- getMCI(cluster, test, fun = 'cor')
#' names(membersL)
#' [1] "members" "MCI"     "sd"      "PCC"     "PCCo"  
#'
#' @import psych
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}

getMCI <- function (groups,  countsL,  adjust.size = FALSE,  
                    fun = c("cor", "BioTIP"), df = NULL   
) 
{
  fun <- match.arg(fun)
  PCC_gene.target = 'zero'
  if (all(is.na(groups))) {
    warning("no loci in any of the state in the list given,  
            please rerun getCluster_methods with a larger cutoff 
            or provide a list of loci")
  } else {
    if (all(sapply(groups,  class) ==  "communities")) {
      membersL = lapply(groups,  membership)
    }
    else if (any(is.na(groups)) & any(sapply(groups,  class) ==  
                                      "communities")) {
      removed = groups[is.na(groups)]
      groups = groups[!is.na(groups)]
      membersL = lapply(groups,  membership)
    } else {
      membersL = groups
    }
    CIl = PCCol = PCCl = sdl = list()
    names(membersL) = names(groups)
    if (is.null(names(groups))) 
      warning("No names provided for \"groups\"")
    if (is.null(names(countsL))) 
      warning("No names provided for \"countsL\"")
    
    loop = names(membersL)
    for (i in loop) {
      test = membersL[[i]]
      if (all(is.na(test))) {
        CI = sdL = PCC = PCCo = NA
      } else {
        test.counts = countsL[[i]]  # gene x sample matrix
        m = lapply(1:max(test),  
                   function(x) subset(test.counts, 
                                      row.names(test.counts) %in% names(test[test ==  x])))
        comple = lapply(1:max(test),  
                        function(x) subset(test.counts,  
                                           !row.names(test.counts) %in% names(test[test ==  x])))
        names(m) = names(comple) = 1:max(test)  
        
        if(fun ==  "cor") {
          PCCo = lapply(names(comple),  function(x) abs(cor(t(comple[[x]]),    
                                                            t(m[[x]]))))          
          PCCo_avg = sapply(PCCo,  function(x) mean(x,  na.rm = TRUE))        
          PCC = lapply(m,  function(x) abs(cor(t(x))))
          PCC_avg = sapply(PCC,  
                           function(x) (sum(x,  na.rm = TRUE) - 
                                          nrow(x))/(nrow(x)^2 - nrow(x)))
        }
        if(fun ==  "BioTIP") {
          if(is.null(df)) {
            # warning("MCI with the 'BioTIP'method is performing local estimation") # suppressed on 01/15/2022
            PCCo_avg = lapply(names(comple),  
                              function(x) avg.cor.shrink(comple[[x]], Y = m[[x]], 
                                                         abs = TRUE, 
                                                         MARGIN = 1,  
                                                         target = PCC_gene.target))
            PCCo_avg = unlist(PCCo_avg) 
            PCC_avg = lapply(m,  function(x) avg.cor.shrink(x, Y = NULL, 
                                                            abs = TRUE, 
                                                            MARGIN = 1,  
                                                            target = PCC_gene.target))
            PCC_avg = unlist(PCC_avg) 
          } else {
            M <- cor.shrink(df,  Y = NULL, 
                            MARGIN = 1,  
                            target = PCC_gene.target)
            PCCo_avg <- array(dim = length(m))
            names(PCCo_avg) <- names(m)
            for(j in 1:length(m)){
              PCCo_avg[j] <- mean(abs(M[rownames(comple[[j]]),  rownames(m[[j]])]))
            }
            PCC_avg <- array(dim = length(m))
            names(PCC_avg)  <- names(m)
            
            for(j in 1:length(m)){
              tmp <- M[rownames(m[[j]]),  rownames(m[[j]])]
              U <- upper.tri(tmp,  diag = FALSE)
              PCC_avg[j] <- mean(abs(U))
            }
          }
        }  
        
        sdL = lapply(m,  function(x) apply(x,  1,  sd))
        if (adjust.size) {
          MCI = mapply(function(x,  y,  z,  w) mean(x) * 
                         (y/z) * sqrt(nrow(w)),  sdL,  PCC_avg,  PCCo_avg,  
                       m)
        } else {
          MCI = mapply(function(x,  y,  z) mean(x) * (y/z),  
                       sdL,  PCC_avg,  PCCo_avg)
        }
      }
      CIl[[i]] = MCI
      sdl[[i]] = sdL
      PCCl[[i]] = PCC_avg
      PCCol[[i]] = PCCo_avg
    }
    names(CIl) = names(sdl) = names(PCCl) = names(PCCol) = names(membersL)
    return(list(members = membersL,  MCI = CIl,  sd = sdl,  
                PCC = PCCl,  PCCo = PCCol))
  }
}


#' @title Get MCI Scores for randomly selected genes
#'
#' @description This function gets the MCI scores for randomly selected features (e.g. transcript ids), 
#'
#' @param len An integer that is the length of genes in the CTS (critical transition signal).
#' 
#' @param samplesL A list of vectors,  whose length is the number of states. Each vector gives the sample names in a state. 
#' Note that the vector s (sample names) has to be among the column names of the R object 'df'.
#' 
#' @param df A numeric matrix or dataframe of numerics, factor or character. 
#' The rows and columns represent unique transcript IDs (geneID) and sample names,  respectively
#' 
#' @param adjust.size A boolean value indicating if MCI score should be adjusted by module size (the number of transcripts 
#' in the module) or not. Default FALSE.
#' 
#' @param B An integer, setting the permutation with \code{B} runs. Default is 1000.
#' 
#' @param fun A character chosen between ("cor", "BioTIP"), indicating where an adjusted 
#'   correlation matrix will be used to calculate the MCI score.   
#' 
#' @param M a pre-calculated global shrunk matrix, can save calculation if working on the same data
#'   for multiple CTS evaluations.
#'
#' @return A numeric matrix indicating the MCI scores of permutation. 
#' The dimension (row X column) of this matrix is the length of \code{samplesL} * \code{B}.
#' 
#' @export
#' 
#' @examples
#' counts = matrix(sample(1:100, 18), 3, 9)
#' colnames(counts) = 1:9
#' row.names(counts) = c('loci1', 'loci2', 'loci3')
#' cli = cbind(1:9, rep(c('state1', 'state2', 'state3'), each = 3))
#' colnames(cli) = c('samples', 'group')
#' samplesL <- split(cli[, 1], f = cli[, 'group'])
#' simMCI = simulationMCI(2, samplesL, counts, B = 2)
#' simMCI
#' #            [,1]      [,2]
#' #state1  2.924194  2.924194
#' #state2 20.877138 20.877138
#' #state3  2.924194  2.924194
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}

simulationMCI  <- function (len,  samplesL,  df,  adjust.size = FALSE,  B = 1000, 
                            fun = c("cor", "BioTIP"), M = NULL ) 
{
  fun <- match.arg(fun)
  PCC_gene.target = 'zero'
  if (is.null(names(samplesL))) 
    stop("please provide names for list countsL")
  countsL = lapply(samplesL,  function(x) df[,  as.character(x)])
  if (is.null(names(countsL))) 
    names(countsL) = names(samplesL)
  # create progress bar
  pb <- txtProgressBar(min = 0,  max = B,  style = 3)
  if (fun ==  "BioTIP") {
    if(is.null(M)) M <- cor.shrink(df, Y = NULL, MARGIN = 1, shrink = TRUE, 
                                   target = PCC_gene.target)
  }
  else M = NULL
  m <- matrix(nrow = length(samplesL), ncol = B)
  for (i in 1:B) {
    setTxtProgressBar(pb, i)
    m[, i] <- getMCI_inner(len, countsL, adjust.size, fun = fun, 
                           PCC_gene.target = PCC_gene.target, M = M)
    Sys.sleep(0.01)
    if (i ==  B) 
      cat("Done!\n")
  }
  row.names(m) = names(countsL)
  return(m)
}


#' @title Calculating MCI Score for randomly selected
#'
#' @description This function calculates random MCI score,  allowing an estimation 
#' of correlation matrix using the Schafer-Strimmer Method for the PCC_in component in the MCI score.
#' 
#' @param members An integer that is the length of genes in the CTS (critical transition signal).
#' 
#' @param countsL A list of subset of the data matrix.  The list length is the number of states.    
#' Each subset of matrix gives the genes in rows and samples in columns.  
#' 
#' @inheritParams  simulationMCI
#' 
#' @param PCC_gene.target A character 'zero' indicating that beetween-gene correlation matrix will be shrunk.
#' towards zero, used only for fun = 'BioTIP'.
#' 
#' @param  M is the overall shrunk correlation matrix, used only for fun = 'BioTIP'. 
#' 
#' @return A vector recording one MCI score per state.
#' 
#' @examples
#' 
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}

getMCI_inner = function(members, countsL,  adjust.size, 
                        fun = c("cor", "BioTIP"),   
                        PCC_gene.target = 'zero',  M = NULL
) 
{
  fun <- match.arg(fun)
  
  random_id = sample(1:nrow(countsL[[1]]), members)
  randomL = lapply(names(countsL),  
                   function(x) countsL[[x]][random_id, ])
  comple = lapply(names(countsL),  
                  function(x) subset(countsL[[x]], 
                                     !row.names(countsL[[x]]) %in% row.names(randomL[[x]])))
  names(randomL) = names(comple) = names(countsL)
  if(fun ==  "BioTIP") {
    PCCo_avg = array(dim = length(countsL))
    names(PCCo_avg) = names(countsL)
    for(i in 1:length(PCCo_avg)) {
      PCCo_avg[i] <- mean(abs(M[rownames(comple[[i]]),  rownames(randomL[[i]])]))
    }
    PCC_avg = array(dim = length(countsL))
    names(PCC_avg) = names(countsL)
    for(i in 1:length(PCC_avg)) {
      PCC_avg[i] = mean(abs(M[rownames(randomL[[i]]),  rownames(randomL[[i]])]))
    }
  } else if(fun ==  "cor") {
    PCCo = lapply(names(comple),  
                  function(x) abs(cor(t(comple[[x]]), t(randomL[[x]]))))
    PCCo_avg = sapply(PCCo, function(x) mean(x, na.rm = TRUE))
    PCC = lapply(randomL, function(x) abs(cor(t(x))))
    PCC_avg = sapply(PCC, 
                     function(x) (sum(x, na.rm = TRUE)-nrow(x))/(nrow(x)^2-nrow(x)))
  }
  sdL = lapply(randomL,  function(x) apply(x, 1, sd))
  
  if(adjust.size){
    MCI = mapply(function(x, y, z, w) mean(x)*(y/z)*sqrt(members),  sdL, PCC_avg, PCCo_avg, members)
  }else{
    MCI = mapply(function(x, y, z) mean(x)*(y/z),  sdL, PCC_avg, PCCo_avg)
  }
  return(MCI)
}



#' @title Plot observed and simulated MCI Scores
#'
#' @description Box plots of observed (red) and simulated MCI scores by boostrapping genes \code{B} times, 
#' with three horizontal lines: the min, max and 2*(max-min) value of the state of interests, or of all values.
#'
#' @param MCI A named vector of max CI scores per state, can be obtained from function \code{\link{getMaxStats}}.
#' 
#' @param simulation A matrix state * number of simulated times,  can be obtained from function \code{\link{simulationMCI}}.
#' 
#' @param las Numeric in {0, 1, 2, 3}; the style of axis labels. Default is 0, meaning labels are parallel. 
#' (link to http://127.0.0.1:21580/library/graphics/html/par.html)
#' 
#' @param order A vector of state names in the customized order to be plotted, set to NULL by default.
#' 
#' @param ylim An integer vector of length 2. Default is NULL.
#' 
#' @param main A character vector. The title of the plot. Default is NULL.
#' 
#' @param which2point A character (or integer) which state's values were used to set up the three horizontal lines. 
#' by default is NULL,  indicating the values of all states will be used.
#' 
#' @export
#' @return Return a box plot of MCI(red) and simulated MCI(grey) scores per state.
#' 
#' @examples
#' MCI = c(1:3); names(MCI) = c('a', 'b', 'c')
#' simMCI = matrix(sample(1:100, 9), 3, 3)
#' row.names(simMCI) = names(MCI)
#' plot_MCI_Simulation(MCI, simMCI)
#'
#' @author Zhezhen Wang \email{zhezhen@@uchicago.edu}

plot_MCI_Simulation <- function(MCI,  simulation,  las = 0, 
                                order = NULL,  ylim = NULL,  main = NULL, 
                                which2point = NULL,  ...)
{
  if(is.null(names(MCI))) stop('make sure elements in "MCI" have names')
  if(!is.null(order)){
    if(any(!order %in% row.names(simulation))) 
      stop('make sure "simulation" has row.names which are in "order"')
    if(any(!row.names(simulation) %in% order)) 
      warning('not every state in "simulation" is plotted,  make sure "order" is complete')
    simulation = simulation[order, ]
  }
  maxpt = max(simulation, MCI, na.rm = TRUE)
  tmp = c(min(simulation, MCI, na.rm = TRUE), maxpt)
  if(is.null(ylim)){
    if(min(simulation, na.rm = TRUE)<maxpt){
      ylim = tmp
    }else{
      ylim = rev(tmp)
    }
  }
  boxplot(t(simulation), col = 'grey', ylab = 'DNB score', 
          axes = FALSE, ylim = ylim, main = main,  pch = 20,  ...)
  
  x = which.max(MCI)
  maxCI = MCI[x]
  
  if(!is.null(order)){
    if(is.null(names(MCI))) stop('make sure "MCI" is named using names in "order"')
  }
  
  axis(2)
  # customize x-axis
  if(is.null(order)){
    stages = row.names(simulation)
  }else{
    stages = order
  }
  x = which(stages ==  names(x))
  axis(side = 1, at = 1:nrow(simulation), labels = stages, las = las)
  points(x, maxCI, col = 'red', pch = 16)
  
  if(is.null(which2point)) {
    abline(h = min(simulation),  col = "grey",  lty = 3)
    abline(h = max(simulation),  col = "grey",  lty = 3)
    abline(h = min(simulation)+ 2*(max(simulation)-min(simulation)),  
           col = "grey",  lty = 2)
  }  else {
    abline(h = min(simulation[which2point, ]),  col = "grey",  lty = 3)
    abline(h = max(simulation[which2point, ]),  col = "grey",  lty = 3)
    abline(h = min(simulation[which2point, ])+ 
             2*(max(simulation[which2point, ])-min(simulation[which2point, ])),  
           col = "grey",  lty = 2)
  }
  
}

#######################################################################
################ new functions,  02/14/2020  ###########################
#######################################################################


#' @title Index of criticality Scoring System with estimated correlation, an updated Ic-score 
#'
#' @description This function calculates the BioTIP score on a given
#' data matrix X (or two matrixes X and Y). It can also calculate the \eqn{I_c} score, if desired.
#' 
#' This approach uses the method outlined by Schafer and Strimmer in
#' "A Shrinkage Approach to Large-Scale Covariance Matrix Estimation
#' and Implications for Functional Genomics" (2005)
#' 
#' This approach is modified to ignore missing values, analogous to how
#'\code{cor(X,  use = "pairwise.complete.obs")} works.
#'
#' The gene-gene correlations are shrunk towards 0,  whereas the
#' sample-sample correlations are shrunk towards their empirical average.
#'
#' @param X A G x S matrix of counts. Rows correspond to genes, 
#' columns correspond to samples.
#'
#' @param method A flag specifying whether to calculate the BioTIP score
#' or the \eqn{I_c} score
#' 
#' @param PCC_sample.target A numeric choose between [0,1] or 
#' a character choose among ('none, 'average',  'zero', 'half'),  
#' indicating whether to turn off shrinkage, or to shrink PCC towards their empirical common average,  
#' zero or 0.5 (for sample-sample correlations).
#' 
#' @param output A string. Please select from 'IndexScore',  'PCCg',  or 'PCCs'. Uses 'IndexScore' by default.
#' 'PCCg' is the PCC between genes (numerator) and 'PCCs' is PCC between samples (denominator).
#'
#' @return A value containing the shrunk BioTIP or non-shrunk \eqn{I_c} score
#'
#'
#' @references
#'
#' @export
#' 
#' @examples
#' ## Generating a data X as coming from a multivariate normal distribution 
#' ## with 10 highly correlated variables, roughly simulating correlated genes.
#' M = matrix(.9, nrow = 10, ncol = 10)
#' diag(M) = 1
#' mu = rnorm(10)
#' X = MASS::mvrnorm(1000, mu, M)
#' dim(X)  #1000 10  
#' 
#' ## Calculating pairwise correlation between 1000 genes; then the mean value
#' ## in two ways, respectively
#' cor_tX = cor(t(X))
#' mean(abs(cor_tX[upper.tri(cor_tX, diag = FALSE)])) # 0.9150228
#' 
#' getIc.new(X, method = "Ic", output = 'PCCg') # 0.9150228
#' getIc.new(X, method = "BioTIP", output = 'PCCg') # 0.8287838
#'
#' ## Using the Index of critical scoring system, in two ways, respectively 
#' (newscore = getIc.new(X, method = "BioTIP"))
#' (oldscore = getIc.new(X, method = "Ic"))
#'
#' @author Andrew Goldstein \email{andrewgoldstein@@uchicago.edu}

getIc.new = function(X,  method = c("BioTIP",  "Ic"),  
                     PCC_sample.target = 1,  ## 12/02/2020
                     output = c('Ic', 'PCCg', 'PCCs')) 
{
  PCC_gene.target = 'zero'
  
  # whether to calculate BioTIP or Ic
  method = match.arg(method)
  # begin ## 12/02/2020
  # PCC_sample.target = match.arg(PCC_sample.target)
  if (class(PCC_sample.target) ==  'numeric') if((PCC_sample.target < 0) | (PCC_sample.target > 1)) 
  { 
    stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
  } else if(class(PCC_sample.target) ==  'character') if(!PCC_sample.target %in% c('none', 'zero',  'average', 'half')) {
    stop("Argument `PCC_sample.target` must be a value between 0 and 1, or a choice of 'none', 'zero',  'average', 'half'")
  }  
  # end ## 12/02/2020
  output <- match.arg(output)
  
  # if using BioTIP,  set shrink = TRUE,  else FALSE
  shrink = (method ==  "BioTIP")
  
  # get numerator and denominator for score
  numerator = avg.cor.shrink(X,  MARGIN = 1,  shrink = shrink,  abs = TRUE,  target = PCC_gene.target)
  denominator = avg.cor.shrink(X,  MARGIN = 2,  shrink = shrink,  abs = FALSE,  target = PCC_sample.target)
  
  if(output ==  'Ic') return(numerator / denominator)
  if(output ==  'PCCg') return(numerator)
  if(output ==  'PCCs') return(denominator)
  
}


#' @title Estimation of average values of correlation 
#'
#' @description This function takes in one (or two) matrix X 
#' (rows are genes, columns are samples) (or Y).
#' It then calculates the average pairwise correlation between genes or samples.
#' This function uses the method outlined by Schafer and Strimmer in
#' "A Shrinkage Approach to Large-Scale Covariance Matrix Estimation
#' and Implications for Functional Genomics" (2005)
#' 
#' @inheritParams cor.shrink
#' 
#' @param abs A flag specifying whether to take the absolute value 
#' before taking the average (used for gene-gene correlations,  
#' not sample-sample correlations) 
#'
#' @return The average pairwise correlation between genes or samples.
#'
#'
#' @references Schafer and Strimmer (2005) "A Shrinkage Approach to Large-Scale 
#' Covariance Matrix Estimation and Implications for Functional Genomics"
#'
#' @export
#' 
#' @examples
#' ## Generating a data X as coming from a multivariate normal distribution
#' ## with 10 highly correlated variables, roughly simulating correlated genes.
#' set.sedd(2020)
#' M = matrix(.9, nrow = 10, ncol = 10)
#' diag(M) = 1
#' mu = rnorm(10)
#' X = t(MASS::mvrnorm(500, mu, M))
#' dim(X)  #10  500, simulating a matrix of 10 genes in rows and 500 samples in columns
#' 
#' ## Calculating pairwise correlation among 10 correlated genes
#' cor_tX = cor(t(X))
#' mean(abs(cor_tX[upper.tri(cor_tX, diag = FALSE)])) # 0.9101072
#' 
#' ## Calculating estimated pairwise correlation among these 10 correlated genes
#' avg.cor.shrink(X, MARGIN = 1,shrink = TRUE, targe = 'zero') # 0.9053253
#' 
#' ## Calculating estimated pairwise correlation
#' ## after adding additional 30 random genes
#' Y = rbind(X, matrix(rnorm(30*500), ncol = 500))
#' dim(Y)  #40  500
#' avg.cor.shrink(Y, MARGIN = 1,shrink = TRUE, targe = 'zero') # 0.04758553
#'
#' Compared with the empirical pairwise correlation 
#' cor_tY = cor(t(Y))
#' mean(abs(cor_tY[upper.tri(cor_tY, diag = FALSE)])) # 0.08547593
#'
#' @author Andrew Goldstein \email{andrewgoldstein@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}
# 
avg.cor.shrink = function(X,  Y = NULL,  MARGIN = c(1,  2),  shrink = TRUE,  
                          abs = FALSE,  target = 0 )  #c('zero',  'average', 'half','none')) ## 12/02/2020
{
  if(target ==  'none')  shrink = FALSE ## begin 12/02/2020
  if (class(target) ==  'numeric') if((target < 0) | (target > 1)) 
    { # make sure we have a valid target for off-diagonal correlations
      stop("Argument `target` must be a value between 0 and 1, or a choice of 'zero',  'average', 'half', 'none'")
    } else if(class(target) ==  'character') if(!target %in% c('none', 'zero',  'average', 'half')) {
      stop("Argument `target` must be a value between 0 and 1, or a choice of 'zero',  'average', 'half','none'")
    }  
  if(MARGIN != 1 & MARGIN != 2) stop("MARGIN must be a choice of 1 or 2.")
  #  MARGIN = match.arg(MARGIN)
  #  target = match.arg(target)
  ## end 12/02/2020
  
  X_cor_shrink = cor.shrink(X = X,  Y = Y,  MARGIN = MARGIN,  shrink = shrink,  target = target)
  if(is.null(Y)){
    U <- upper.tri(X_cor_shrink,  diag = FALSE)
    X_cor_shrink <- X_cor_shrink[U] 
  }
  if (abs ==  TRUE) {
    res = mean(abs(X_cor_shrink), na.rm = TRUE)
  } else {
    res = mean(X_cor_shrink, na.rm = TRUE)
  }
  
  return(res)
}


#'
#' @description This function is the "work-horse" of the BioTIP function. 
#' It takes in a matrix X (rows are genes, columns are samples), and an optional matrix Y of the same form. 
#  It then calculates the average pairwise correlation between genes or samples. 
#  There are options to choose what constant value to shrink the off-diagonal correlations to. 
#  This target must be a value between 0 and 1. 
#' It then calculates the average pairwise correlation between genes or samples. 
#' This method adopts the method outlined by Schafer and Strimmer (2005). 
#' 
#' @param X A G1 x S matrix of counts. Rows correspond to genes, 
#' columns correspond to samples.
#'
#' @param Y A G2 x S matrix of counts. Rows correspond to genes, 
#' columns correspond to samples. By default is NULL.
#'
#' @param MARGIN An integer indicateing whether the rows (1,  genes) 
#' or the columns (2,  samples) to be calculated for pairwise correlation.
#'
#' @param shrink A flag specifying whether to shrink the correlation or not.
#' If the parameter 'target' is 'none', shrink will be set to 'FALSE.'
#'
#' @param target A number choose between 0 and 1, or a character among ('zero', 'average',  'half', 'none'),  
#' indicating whether to turn off 'shrink', or to shrink towards zero (for gene-gene correlations),   
#' shrink towards their empirical common average,  
#' one or 0.5 (for sample-sample correlations).
#'
#' @return The pairwise correlation between genes or samples.
#' If Y ==  NULL, a G1 x G1 matrix is returned; otherwise, a G1 x G2 matrix is returned.
#'
#' @references Schafer and Strimmer (2005) "A Shrinkage Approach to Large-Scale 
#' Covariance Matrix Estimation and Implications for Functional Genomics"
#'
#' @export
#'
#' @examples
#' require(MASS)
#' ## Generating a data X as coming from a multivariate normal distribution 
#' ## with 10 highly correlated variables, roughly simulating correlated genes.
#' set.seed(2020)
#' M = matrix(.9, nrow = 10, ncol = 10)
#' diag(M) = 1
#' mu = rnorm(10)
#' X = t(MASS::mvrnorm(500, mu, M))
#' dim(X)  #10 500
#' 
#' ## Calculating pairwise correlation among 10 genes
#' cor_tX = cor(t(X))
#' mean(abs(cor_tX[upper.tri(cor_tX, diag = FALSE)])) # 0.9109133
#' 
#' ## Calculating estimated pairwise correlation among 500 samples
#' cor.matrix <- cor.shrink(X, MARGIN = 1, shrink = TRUE, targe = 0) 
#' dim(cor.matrix)   #[1] 10  10
#' mean(upper.tri(cor.matrix, diag = FALSE))  # 0.45
#' 
#' ## Calculating estimated pairwise correlation again, 
#' ## after adding additional 100 random genes
#' Y = rbind(X, matrix(rnorm(500*100), ncol = 500))
#' dim(Y)  #110 500
#' cor.matrix <- cor.shrink(Y, MARGIN = 1, shrink = TRUE, targe = 0.5) 
#' dim(cor.matrix)   #[1] 110  110
#' mean(upper.tri(cor.matrix, diag = FALSE))  # 0.4954545
#' 
#' ## alternatively, you can run
#' Y = matrix(rnorm(500*100), ncol = 500)
#' cor.matrix <- cor.shrink(X, Y, MARGIN = 1, shrink = TRUE, targe = 0.5) 
#' dim(cor.matrix)   #[1] 120  120
#' mean(upper.tri(cor.matrix, diag = FALSE))  # 0.4954545
#' 
#' @author Andrew Goldstein \email{andrewgoldstein@@uchicago.edu}; Xinan H Yang \email{xyang2@@uchicago.edu}
# 

cor.shrink = function(X,  Y = NULL,  MARGIN = c(1,  2),  shrink = TRUE,  
                      target = 0) ## 12/02/2020
{
  # We have the two following reasons to shut off the parameter target='average' here:
  # Theoretically, the 'TARGET D' outlined by Schafer and Strimmer (2005) can't be fed by the estimated average.
  # Practically, the shrinkage towards average generates an estimated matrix, 
  # whose average value remains the same as that of its observation matrix.
  if(target ==  'average')  shrink = FALSE ## 12/02/2020 
  if(target ==  'none')  shrink = FALSE ## begin 12/02/2020
  # MARGIN = match.arg(MARGIN)
  
  # target = match.arg(target) 
  if (class(target) ==  'numeric') if((target < 0) | (target > 1)) { # make sure we have a valid target for off-diagonal correlations
    stop("Argument `target` must be a value between 0 and 1, or a choice of 'zero',  'average', 'half', 'none'")
  } else if(class(target) ==  'character') if(!target %in% c('zero',  'average', 'half', 'none')) {
    stop("Argument `target` must be a value between 0 and 1, or a choice of 'zero',  'average', 'half', 'none")
  }
  ## end 12/02/2020
  
  dim_X = dim(X) 
  dim_Y = dim(Y) 
  
  
  Y.exist = FALSE
  
  if(!is.null(Y)){
    if(MARGIN ==  1) {
      X <- rbind(X, Y)
      Y.exist = TRUE  # a flag for outputs
    } else {
      X <- cbind(X, Y)
    }
  } 
  
  # center and scale X by mean and SD (using wither column or row means/SDs)
  X_means = apply(X,  MARGIN = MARGIN,  mean,  na.rm = TRUE)
  X_sds = apply(X,  MARGIN = MARGIN,  sd,  na.rm = TRUE)
  X_sds[X_sds ==  0] = 1 # in case where variable doesn't vary
 
  X_std = sweep(sweep(X,  MARGIN = MARGIN,  STATS = X_means,  FUN = '-'), 
                MARGIN = MARGIN,  STATS = X_sds,  FUN = '/')
  
  # From now on Y is a new internal matrix regardless the input Y ### 
  rm(Y) #update on 9/29/2020
  
  Y = !is.na(X_std) # contains pattern of non-missing values
  X_std[!Y] = 0 # set missing to 0 to not add to (t)crossprod (but # of not missing encoded in Y)
  
  # calculate quantities needed
  if (MARGIN ==  1) {
    XtX = tcrossprod(X_std)
    X2tX2 = tcrossprod(X_std^2)
    YtY = tcrossprod(Y)
  } else {
    XtX = crossprod(X_std)
    X2tX2 = crossprod(X_std^2)
    YtY = crossprod(Y)
  }
  
  # calculate empirical correlation
  X_cor = XtX / (pmax(1,  YtY - 1, na.rm = TRUE))
  X_cor_shrink = X_cor # initialize with empirical correlation
  
  U = upper.tri(XtX,  diag = FALSE) # index on U to get upper triangular part of matrix
  
  if (shrink) { # if we're shrinking,  find lambda and shrink towards specified target
    numerator = sum(((YtY[U] * X2tX2[U]) - (XtX[U])^2) / ((pmax(1,  YtY[U] - 1))^3), na.rm = TRUE)
    
    # BEGIN translate character to numeber
    if(class(target) ==  'character'){ #12/02/2020
      if (target ==  'zero') {  
        target = 0
      } else {
        if(target ==  'half') {
          target = 0.5
        } else {
          target = mean(X_cor[U])
        }
      }
    }    #12/02/2020  
    denominator = sum((X_cor[U]-target)^2)   #update on 9/29/2020
    
    
    lambda = ifelse(shrink,  max(0,  min(1,  numerator / denominator), na.rm = TRUE),  0)
    target_cor = matrix(target, nrow = nrow(X_cor), ncol = ncol(X_cor))
    
    diag(target_cor) = 1
    
    X_cor_shrink = (lambda * target_cor) + ((1 - lambda) * X_cor)
  }
  #   cat(dim(X_cor), '\t')
  
  # if(Y.exist){ ## 12/02/2020 removed this condition to report consistent output
  #   X_cor_shrink = X_cor_shrink[1:dim_X[MARGIN], (dim_X[MARGIN] + 1):(dim_X[MARGIN] + dim_Y[MARGIN])]  #update on 9/29/2020
  # }   
  
  return(X_cor_shrink)
}



#' @title Density plot the leading two distance between any two states from random scores of all states in a system.
#'
#' @description Generate a density plot of Ic score (orBioTIP score) from a simulation, 
#' which is the distance between the first-larget and the second-largest random scores. 
#' This is an alternative method to estimate the significance of an observed BioTIP (or Ic) score in a system. 
#' This measurement makes more sense to evaluate random scores of sample-label shuffling, 
#' in which the nature sample-sample correlation within a phenotypic state (or cell subpopulation) was removed. 
#'
#' @inheritParams plot_Ic_Simulation 
#' 
#' @param xlim An integer vector of length 2. Default is NULL.
#' 
#' @param which2point A character (or integer) which state's values were used to set up the three horizontal lines. 
#' by default is NULL, indicating the values of all states will be used.
#' 
#' @export
#' 
#' @return Return a P value and a plot of the observed Ic (red) and simulated Ic (grey) scores per state.
#' 
#' @author Xinan H Yang \email{xyang2@@uchicago.edu}
#' 
#' @examples
#' sim = matrix(sample(1:10, 9), 3, 3)
#' row.names(sim) = paste0('state', 1:3)
#' Ic = c('state1' = 3.4, 'state2' = 15.6, 'state3' = 2)
#' plot_SS_Simulation(Ic, sim)
#' 


plot_SS_Simulation <- function(Ic,  simulation,  las = 0,  
                               xlim = NULL,  ylim = NULL,  order = NULL,  
                               main = "1st max - 2nd max",     
                               ylab = "Density")
{
  if (any(is.null(names(Ic)))) 
    stop("Please provide name for vector \"Ic\" ")
  if (!identical(names(Ic),  row.names(simulation))) 
    Ic = Ic[match(row.names(simulation),  names(Ic))]
  if (!is.null(order)) {
    if (any(!names(Ic) %in% order)) 
      warning("not all states in Ic is plotted")
    if (any(!order %in% names(Ic))) 
      stop("make sure \"Ic\" is named using names in \"order\"")
    toplot = toplot[order,  ]
  }
  # if(!is.null(which2point)) {
  #   if (!class(which2point) %in% c("integer", "character")) 
  #     stop("'which2poin' should be an integer or characer.")
  #   if (class(which2point) ==  "integer" & which2point > length(Ic)) 
  #     stop (paste("'which2point' should be small than",  length(Ic)))
  #   if (class(which2point) ==  "character" & !which2point %in% names(Ic)) 
  #     stop (paste("'which2point' should one of",  names(Ic)))
  # }
  
  diff_Ic <- apply(simulation,  MARGIN = 2,   
                   function(x) sort(x,  decreasing = TRUE)[2:1])
  diff_Ic <- apply(diff_Ic,  MARGIN = 2,  diff)
  density_diff = density(diff_Ic)  # # Kernel Density
  if(is.null(ylim))  # correct 09/29/2020
    ylim = c(0,  .1 + max(density_diff$y))
  if(is.null(xlim)) {  # correct 09/29/2020
    xlim = c(-.05,  .05) + c(min(density_diff$x),  max(density_diff$x))
  }
  plot(density_diff,  type = 'l',  lwd = 2,  col = 'black',  
       main = main,  ylab = paste(ylab, "Density"), 
       xlim = xlim,  
       ylim = ylim,  
       cex.main = 1.2,  cex.lab = 1.2)
  
  v <- diff(sort(Ic,  decreasing = TRUE)[2:1])
  
  abline(v = v,  col = 'red',  lwd = 2,  lty = 2)
  P.value <- round(mean(diff_Ic >= v), 3)
  legend("topright",  
         legend = paste0("p = ",   
                         P.value 
         ), 
         text.col = "red",  bty = 'n',  cex = 1.5)
  
  return(P.value)  ## add by xy on 6/23/2020
}


#' @title Identifying the 'Biomodule', the candidate critical transition signal transcripts
#'
#' @description This function reports the 'biomodule', which is the module with
#'   the maximum Module Critical Index (MCI) scores for each state. Each state
#'   can have multiple modules (groups of subnetworks derived from the function
#'   \code{\link{getCluster_methods}}). This function runs over all states.
#'
#' @param membersL A list of vectors with unique sample ids as cluster names. The length
#'   of this list is equal to the number of states in the study. This can be the
#'   first element of the output from function \code{getMCI} or the output from
#'   \code{getCluster_methods}, see Examples for more detail.
#'   
#' @param idL A list of numeric vectors with unique cluster numbers as names.
#'   This can be the first element of the output from function \code{getMaxMCImember}.
#'   
#' @param whoisnext A vector of who is the next state
#' in the list given by the \code{memebersL} to exract.
#'   
#' @param which.next An integer indicating how many modules with top MCI scores to output. 
#' By default is 2 meaning the 2nd id.
#'
#' @return A nested list whose length is the length of the input object
#'   \code{membersL}.  Each internal list contains at least two objects: one object is
#'   the vector of the reported biomodule per state, and the other object is a list of
#'   transcript IDs (each defines the biomodule of maximum score per state) across states.
#'   When n>1, there will be additional object(s) to record the list of transcrit IDs 
#'   (each defines the biomodule of the 2nd , the 3rd , ..., maximum score per state)
#' @export
#' A vector of MCI score for each state of interest
#' 
#' @author Xinan H Yang \email{xyang2@@uchicago.edu}
#' 
#' @examples
#' set.seed(2020)
#' test = list('state1' = matrix(rnorm(50,0,1), 10, 5), 
#'             'state2' = matrix(rnorm(30,0,3), 10, 3), 
#'             'state3' = matrix(rnorm(40,0,1), 10, 4))
#' 
#' ## Assign colnames and rownames to the matrix
#' for(i in 1:length(test)){
#'   colnames(test[[i]]) = paste0('Cell_',i,1:ncol(test[[i]]))
#'   row.names(test[[i]]) = paste0('g',1:10)
#' }
#' 
#' gene.cluster = list(rep(1:2, 5), c(rep(1:3,3),2), rep(1:5,2))
#' names(gene.cluster) = names(test)
#' for(i in names(gene.cluster)){
#'   names(gene.cluster[[i]]) = paste0('g',1:10)
#' }
#' 
#' membersL <- getMCI(gene.cluster, test)
#' names(membersL)
#' # [1] "members" "MCI"     "sd"      "PCC"     "PCCo" 
#' 
#' ## A list of index of interested gene.cluster IDs per state
#' idxL = list(state1 = c(1,2), state1 = integer(0), state3 = c(3,4,1))
#' 
#' ## Extract the score of the 2nd (by default) gene,cluster ID in the 'whoisnext' state
#' getNextMaxStats(membersL[['MCI']], idxL, whoisnext = 'state3')
#' ## Extract the score of the 3rd gene,cluster ID in the 'whoisnext' state
#' getNextMaxStats(membersL[['MCI']], idxL, whoisnext = 3, which.next = 3)
#' 
getNextMaxStats <- function(membersL, idL = maxMCIms[['idx']], whoisnext, which.next = 2)
{
  score <- array(dim = length(whoisnext))
  names(score) <- whoisnext
  idx <- lapply(whoisnext, function(x) idL[[x]][which.next])
  for(i in 1:length(whoisnext))
  {
    score[[i]] <- membersL[[whoisnext[i]]][idx[[i]]]
  }
  
  return(score) 
}

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
#' @author Xinan H Yang \email{xyang2@@uchicago.edu}

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
