#'
#' Classifying Transcriptomes into Biotypes
#'
#' @aliases getBiotypes
#'
#' @description
#' The purpose of the getBiotypes function is to class both coding and noncoding
#' transcripts into biotypes using the most recent GENCODE annotations. This
#' tool can also be used to define potential lncRNAs, given an available genome
#' transcriptome assembly (a gtf file) or any genomic loci of interest.
#'
#' @param full_gr A GRanges object of coding and noncoding transctipts. Unique
#'   identifications for each column must be assigned.More details can be found
#'   in the GRanges package.
#' @param gencode_gr This GRanges object contains a GENCODE reference
#'   annotation.It must have a column of biotypes.
#' @param intron_gr A GRanges object containing the coordinates of introns.For
#'   details see GRanges package.
#' @param minoverlap Detects minimum overlap between two IRanges objects.
#'   Details Overlap arguments are included in the IRanges package.
#'
#' @details For details of findOverlaps, type.partialOverlap, type.50Overlap
#' type.toPlot, queryhits, and subjecthits see
#' [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
#' [IRanges](https://www.bioconductor.org/packages/release/bioc/html/IRanges.html),
#' and [BiocManager](http://bioconductor.org/install/index.html).
#'
#' @return Returns the classified transcriptome biotypes.
#'
#' @source [Refrence GRCh37 genome](https://www.gencodegenes.org/human/release_25lift37.html)
#' for details on gtf format visit [ensemble](https://useast.ensembl.org/info/website/upload/gff.html)
#' @import GenomicRanges
#'
#' @references
#'
#' Wang, Z. Z., J. M. Cunningham and X. H. Yang (2018).'CisPi: a transcriptomic score for disclosing cis-acting disease-associated lincRNAs.'
#' Bioinformatics34(17): 664-670'
#'
#' @examples
#' #Input datasets locally
#' library("GenomicRanges" , "IRanges")
#'
#' data("gencode")
#' data("intron")
#' data("ILEF")
#' gencode_gr = GRanges(gencode)
#' ILEF_gr = GRanges(ILEF)
#' cod_gr = GRanges(cod)
#'
#' getBiotypes(ILEF_gr,gencode_gr)
#'
#' \dontrun{getBiotypes(intron_gr)}
#'
#' @note
#' Replace the path_file when loading data locally to the data directory.
#'
#' @import GenomicRanges IRanges GenomeInfoDb stats
#' @importFrom stats aggregate
#' @export

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

#'
#' Finding Overlaps in Coding Regions
#'
#' @description
#' The getReadthrough function is used to find long transcripts that cover more
#' than two coding regions for gene regions of interst.
#'
#' @param gr A GRanges object that shows the start and end loci on genome.
#' @param cod_gr A GRanges object that contains coding regions. For details
#'   please visit /R/data.R.
#'
#' @details For details of findOverlaps, type.partialOverlap, type.50Overlap
#'   type.toPlot, queryhits, readthrough and subjecthits see,
#'   [GenomicRanges](https://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html),
#'   [IRanges](https://www.bioconductor.org/packages/release/bioc/html/IRanges.html),
#'    and [BiocManager](http://bioconductor.org/install/index.html).
#'
#' @return Returns the classified transcriptome biotypes.
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
#' library("GenomicRanges" , "IRanges")
#'
#' data("gencode")
#' data("ILEF")
#' data("cod")
#' gencode_gr = GRanges(gencode)
#' ILEF_gr = GRanges(ILEF)
#' cod_gr = GRanges(cod)
#'
#' getReadthrough(ILEF_gr,cod_gr)
#'
#' \dontrun{getReadthrough(cod_gr)}
#'
#' @note
#' Replace the path_file when loading data locally to the data directory.
#'
#' @import GenomicRanges IRanges GenomeInfoDb
#' @importFrom stats aggregate
#' @export


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

#' The length of a string (in characters).
#'
#' @param df a count matrix with unique loci row names:ID X column names:samples.
#' @param samplesL a list of characters with stages as names.
#' @param cutoff numeric value, if < 1 automaticlly goes to select top x% transcripts
#'    if > 1 goes to select by x-fold more than the selected method (which is
#'    either the reference, other stages or pervious stage) default 0.01.
#' @param method select from 'reference','other', or 'previous'.
#' for 'reference', the reference has to be the first
#' for 'previous', make sure sampleL is in the right order from benign to malign
#' default uses 'other'
#'
#' @return numeric vector giving number of characters in each element of the
#'   character vector. Missing strings have missing length.
#' @export
#' @examples
#' GSE6136 = read.table("C:/Users/benjk/Desktop/GSE6136_matrix.txt", header = TRUE, comment = '!')
#' row.names(GSE6136) = GSE6136$ID_REF
#' GSE6136 = GSE6136[,-1]
#' cli = read.delim("C:/Users/benjk/Desktop/GSE6136_cli.txt", head = F)
#' cli = t(cli)
#' colnames(cli) = str_split_fixed(cli[1,],'_',2)[,2]
#' cli = cli[-1,]
#' cli = data.frame(cli)
#' cli[,'cell-type:ch1'] = str_split_fixed(cli$characteristics_ch1.1,': ',2)[,2]
#' cli[,'Ig clonality:ch1'] = str_split_fixed(cli$characteristics_ch1.3,': ',2)[,2]
#' colnames(cli)[colnames(cli) == 'cell-type:ch1'] = 'group'
#' cli$Row.names = cli[,1]
#'
#' dat <- GSE6136
#' df <- log2(dat+1)
#' tmp <- names(table(cli$group))
#' sampleL <- split(cli[,1],f = cli$group)
#' test <- sd_selection(df,sampleL,0.01)
#' View(test)
#'
#' @author Zhezhen Wang and Biniam Feleke

sd_selection <- function(df,samplesL,cutoff = 0.01, method = 'other'){
  if(is.null(names(samplesL))) stop('please provide name to samplesL')
  tmp = names(samplesL)
  samplesL = lapply(samplesL,as.character)
  test2 = sapply(tmp, function(x) apply(df[,samplesL[[x]]],1,sd,na.rm = T))

  if(method == 'reference'){
    ref = as.character(sampleL[[1]])
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

#' get an igraph object based on Pearson Correlation Coefficiency(PCC)
#'
#' @description get a igraph object based on Pearson Correlation Coefficiency(PCC)
#'
#' @param optimal a list of count matrix
#' @param p a numeric cutoff
#' @return a list of an igraph object
#' @export
#' @examples
#' igraphL <- getNetwork(test, p=1)
#'
#' @importFrom stringr psych igraph
#'
#' @author Zhezhen Wang and Biniam Feleke

getNetwork <- function(optimal, p = 0.05){
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
    tmp = lapply(1:nrow(test),function(x) test[x,test.p[x,]<p])
    tmp_name = lapply(1:nrow(test),function(x) which(test.p[x,]<p))
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

#' getCluster
#'
#' @param igraphL a list of count matrix or a list of igraph object
#' @param steps number of steps
#' @return A list of objects
#' @export
#' @examples
#' getCluster(igraphL, steps = 4)
#'
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

#' get clusters of nodes by clustering methods
#'
#' @param igraphL a list of count matrix or a list of igraph object
#' @param method a character need to be select from 'rw', 'hcm', 'km', 'pam', or
#'   'natrual'. Uses 'rw' by default for 'hcm', using complete. for 'km' using
#'   'euclidean'.
#' @param cutoff a numeric value, default NULL
#' @param countsL a list of count matrix
#' @return A list of vertex id vectors that can be used in plot function in igraph package as 'mark.groups' parameter
#' @export
#' @examples
#' cluster <- getCluster_methods(igraphL)
#' head(cluster)
#'
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
    groups = lapply(1:length(igraphL), function(x) components(igraphL[[x]])$membership)

  }else(stop('please select from "rw", "hcm","km", "pam", "natrual" as method'))
  return(groups)
}

# step3: calulate CI

#' get CI score
#'
#' @description get CI score
#'
#' @param groups A list of vertex id vectors
#' @param countsL a list of numeric count matrix
#' @param plot a boolean to decide if plot a bar plot of CI scores or not. Default TRUE
#' @param adjust.size a boolean to decide if CI score should be adjust by size or not. Default FALSE
#' @param ylim a vector needed if the output barplots need to be on the same y scale
#' @param nr the number of rows to plot
#' @param nc the number of column to plot
#' @param order the order of the barplot. default NULL which is using the input list order
#' @return a list of CI score and their components
#' @export
#' @examples
#' par(mfrow <- c(1,5))
#' membersL_noweight <- getCI(cluster,test,adjust.size = F,ylim = c(0,8))
#' head(membersL_noweight)
#' maxCIms <- getMaxCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =3)
#' head(maxCIms)
#'
#' @importFrom psych igraph
#' @author Zhezhen Wang
#

getCI <- function(groups,countsL,plot = TRUE, adjust.size = FALSE,ylim = NULL,nr=1,nc = length(countsL),order = NULL){
  if(any(sapply(groups,class) == "communities")){
    membersL = lapply(groups,membership)}
  else{
    membersL = groups
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
        barplot(CI,legend = paste0('#',names(m),'=',sapply(m,nrow)),col = rainbow(length(m), alpha = 0.3),
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

#' get the index and members with the maximum CI score
#'
#' @description get the index and members with the maximum CI score
#'
#' @param membersL a list of characters, the first element of output from function getCI
#' @param CIl alist of numerics, the secned element of output from function getCI
#' @param minsize a numeric value of minimum size allowed for a cluster
#' @return a list of index and members
#' @export
#' @examples
#' membersL_noweight <- getCI(cluster,test,adjust.size = F,ylim = c(0,8))
#' maxCIms <- getMaxCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =3)
#' head(maxCIms)
#'
#' @author Zhezhen Wang

getMaxCImember <- function(membersL,CIl,minsize = 1){
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


