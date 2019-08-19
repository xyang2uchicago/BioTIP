#' A chr21 data from GENCODE GRCh37
#'
#' A reference human GRCh37 genome dataset from GENCODE.
#' \url\{https://www.gencodegenes.org/human/release_25lift37.html}. Included
#' variables are as follows:
#' @format A data frame of chr21 with 2619444 rows and 9 variables:
#' \describe{
#'   \item{Rows}{Rows, include chromosome numbers}
#'   \item{Columns}{Columns include seqname, source, feature, start , end,
#'   score, strand, frame and attribute}
#'   \item{cod_gr}{cod_gr is a subset of sample data 'gencode_gr'}
#' }
#'
"gencode"

#' Coding transcriptome in chr21 dataset
#'
#' A dataset containing chromosomes in the genome regions of interest.
#'  The variables are as follows:
#'
#' @format A data frame with 659327 rows and 5 variables:
#' \describe{
#'   \item{seqnames}{chromosome names (chr1,chrM,chr21)}
#'   \item{ranges}{chromosome ranges on the genome(167684657--167729508)}
#'   \item{strand}{specific strand of the genomic location (+,-,*)}
#'   \item{name}{internal assigned names(uc001aaa.3_intron_0_0_chr1_12228_f)}
#'   \item{score}{score not used in this data set(0)}
#' }
"intron"

#' Chromosome ranges of chr21 dataset
#'
#' A dataset containing chromosomes in the genome regions of interest for 137
#' chromosome. The variables are as follows:
#'
#' @format A data frame of chr21 with 137 rows and 4 variables:
#' \describe{
#'   \item{seqnames}{chromosome names (chr1,chr6,chr21)}
#'   \item{start}{gene read start position (167684657,167729508)}
#'   \item{end}{end of gene read position (15710335,43717938)}
#'   \item{width}{width of gene}
#'   \item{strand}{specific strand of the genomic location (+,-,*}
#'   \item{Row.names}{name of the data rows(A1BG,vawser)}
#' }
"ILEF"

#' cod dataset
#'
#' A subset of gencode_gr extracted as:
#'   cod <- subset(gencode_gr, biotype == 'protein_coding')
#'
#' @format A dataframe with 3 data columns and 1 metadata column.
#' \describe{
#'   \item{seqnames}{chromosome names (chr21)}
#'   \item{ranges}{ranges}{chromosome ranges on the genome (10906201-11029719)}
#'   \item{strand}{specific strand of the genomic location (+,-,*)}
#' }
"cod"


#' GSE6136 matrix dataset
#'
#' A gene annotation samples with their corresponding geneId extracted from GENCODE database.
#' A python script was then run to extcat GSE6136_cli dataset. The script is included in the R/data_raw folder.
#'
#' @format A dataframe with 22690 columns and 27 rows column.
#' \describe{
#'    \item{GSM142398-GSM142423}{Names of gene interest}
#'    \item{ID_REF}{Reference ID of the target genes}
#' }
#'@source \url{https://www.gencodegenes.org/human/}
"GSE6136_matrix"

#' GSE6136 cli dataset
#'
#' A gene annotation samples with their corresponding geneId extracted from GENCODE database.
#' A python script was then run to extcat GSE6136_cli dataset. The script is included in the R/data_raw folder.
#'
#' @format A dataframe with 22690 columns and 27 rows column.
#' \describe{
#'    \item{GSM142398-GSM142423}{Names of gene interest}
#'    \item{Rows}{Summary of GSM genes}
#' }
#'@source \url{https://www.gencodegenes.org/human/}
"GSE6136_cli"
