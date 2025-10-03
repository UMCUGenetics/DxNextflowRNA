#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(FRASER)
  library(biomaRt)
  library(utils)
  library(argparse)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

#' Parse command-line arguments for FRASER run
#'
#' This function defines and parses the command-line arguments required to run
#' FRASER. It uses the \pkg{argparse} parser and returns all arguments as a list.
#'
#' @param input Character vector. Input BAM files.
#' @param refset Character vector. Reference set BAM files.
#' @param prefix Character scalar. Output prefix for results.
#' @param gtf Character scalar. Path to gene annotation file (GTF).
#' @param paired Logical. Whether data are paired-end reads (`TRUE` by default).
#' @param threads Integer. Number of parallel threads to use (default: 12).
#' @param minExpressionInOneSample Integer. Minimum expression in at least one
#'   sample (default: 20).
#' @param quantileMinExpression Integer. Minimum expression based on quantile
#'   (default: 10).
#' @param quantile Numeric. Quantile threshold for expression filtering
#'   (default: 0.95).
#' @param minDeltaPsi Numeric. Minimum delta PSI to consider (default: 0.05).
#' @param filter Logical. If `TRUE`, filter data points rather than masking them
#'   (reduces memory usage). Default: `TRUE`.
#' @param version Character. Prints the FRASER package version and exits.
#'
#' @return A list of parsed arguments as returned by
#'   \code{argparse::ArgumentParser$parse_args()}.
#'
#' @examples
#' \dontrun{
#' args <- get_arguments()
#' }
#'
#' @seealso \code{\link[argparse]{ArgumentParser}}
#' @export
 get_arguments <- function(){
  parser <- ArgumentParser(description="Run FRASER")

  parser$add_argument(
    "-i", "--countTable",
    metavar = "countTable",
    help    = "Input countTable with sample metadata"
  )
  parser$add_argument(
    "--junctionCounts",
    metavar = "junctionCounts",
    help    = "Precalculated junction counts"
  )
  parser$add_argument(
    "--spliceCounts",
    metavar = "spliceCounts",
    nargs   = "+",
    help    = "Precalculated splice site counts"
  )
  parser$add_argument(
    "--prefix",
    metavar = "prefix",
    help    = "Output prefix"
  )
  parser$add_argument(
    "--gtf",
    metavar = "gtf",
    help    = "Gene annotations"
  )
  parser$add_argument(
    "--paired",
    metavar = "paired",
    help    = "Paired-end <true|false>",
    type    = "logical",
    default = TRUE
  )
  parser$add_argument(
    "-p", "--threads",
    metavar = "threads",
    help    = "Number of parallel threads",
    type    = "integer",
    default = 12
  )
  parser$add_argument(
    "--minExpressionInOneSample",
    metavar = "minExpressionInOneSample",
    help    = "minExpressionInOneSample",
    type    = "integer",
    default = 20
    )
  parser$add_argument(
    "--quantileMinExpression",
    metavar = "quantileMinExpression",
    help    = "quantileMinExpression",
    type    = "integer",
    default = 10
    )
  parser$add_argument(
    "--quantile",
    metavar = "quantile",
    help    = "quantile",
    type    = "integer",
    default = 0.95
    )

  parser$add_argument(
    "--minDeltaPsi",
    metavar = "minDeltaPsi",
    help    = "minDeltaPsi",
    type    = "integer",
    default = 0.05
    )
  parser$add_argument(
    "--filter",
    metavar = "filter",
    help    = "filter data points, rather than masking them (reduces memory usage)",
    type    = "logical",
    default = TRUE
    )
  parser$add_argument(
    "-v", "--version",
    action  = "version",
    version = as.character(packageVersion("FRASER"))
  )


  return(parser$parse_args())
}


#' Create a sample table for FRASER input
#'
#' This function generates a sample table from input and reference BAM files,
#' including sample IDs, BAM file paths, and paired-end information.
#'
#' @param sample Character vector. Paths to input BAM files.
#' @param ref Character vector. Paths to reference BAM files.
#' @param paired_end Logical. Whether the data are paired-end reads.
#'   Currently always set to \code{TRUE} in the output.
#'
#' @return A \code{DataFrame} containing:
#' \itemize{
#'   \item \code{sampleID} – BAM file basename without extension
#'   \item \code{bamFile} – full path to BAM file
#'   \item \code{pairedEnd} – logical indicator for paired-end reads
#' }
#'
#' @examples
#' \dontrun{
#' sample_table <- create_sample_table(
#'   sample = c("sample1.bam", "sample2.bam"),
#'   ref = c("ref1.bam", "ref2.bam"),
#'   paired_end = TRUE
#' )
#' }
#'
#' @seealso \code{\link[S4Vectors]{DataFrame}}
#' @export
create_sample_table <- function(sample, ref, paired_end){
  inputs <- c(sample, ref)

  Sampletable <- data.frame(
    sampleID = tools::file_path_sans_ext(basename(inputs)),
    bamFile = inputs,
    stringsAsFactors = FALSE
  )

  Sampletable$pairedEnd <- paired_end



  return(DataFrame(Sampletable))
}


#' Run the FRASER pipeline
#'
#' This function runs the complete FRASER workflow on RNA-seq data:
#' counting, calculating splicing metrics, filtering, annotating with
#' transcript and gene information, and fitting the statistical model.
#'
#' @param sampleTable A \code{DataFrame} containing sample metadata and BAM
#'   file paths, as returned by \code{create_sample_table()}.
#' @param threads Integer. Number of parallel threads to use for computation.
#' @param prefix Character scalar. Prefix for output files (currently unused).
#' @param txdb A \code{TxDb} object with transcript annotations.
#' @param orgdb An organism-specific \code{OrgDb} object for gene annotations.
#' @param minExpressionInOneSample Integer. Minimum read count required in at
#'   least one sample (default: 20).
#' @param quantileMinExpression Integer. Quantile-based expression filter
#'   threshold (default: 10).
#' @param quantile Numeric. Quantile threshold for expression filtering
#'   (default: 0.95).
#' @param minDeltaPsi Numeric. Minimum absolute change in PSI required to keep
#'   a splice event (default: 0.05).
#' @param filter Logical. If \code{TRUE}, data points are filtered rather than
#'   masked (reduces memory usage). Default: \code{TRUE}.
#'
#' @return A \code{FraserDataSet} object containing the processed data and model
#'   results from the FRASER pipeline.
#'
#' @details
#' The pipeline consists of the following steps:
#' \enumerate{
#'   \item Initialize a \code{FraserDataSet} with sample metadata.
#'   \item Count RNA-seq splice junction reads.
#'   \item Calculate PSI values.
#'   \item Filter events by expression and variability criteria.
#'   \item Annotate splice sites with gene and transcript information.
#'   \item Fit the FRASER model to detect aberrant splicing.
#' }
#'
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#'
#' sample_table <- create_sample_table(
#'   sample = c("sample1.bam"),
#'   ref = c("ref1.bam"),
#'   paired_end = TRUE
#' )
#'
#' fds <- run_fraser(
#'   sampleTable = sample_table,
#'   threads = 8,
#'   prefix = "fraser_run",
#'   txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   orgdb = org.Hs.eg.db
#' )
#' }
#'
#' @seealso \code{\link[FRASER]{FraserDataSet}}, \code{\link[FRASER]{FRASER}}
#' @export
run_fraser <- function(fds, threads, prefix, txdb, orgdb, minExpressionInOneSample = 20, quantileMinExpression = 10,
                       quantile = 0.95, minDeltaPsi = 0.05, filter = TRUE){
  bp <- MulticoreParam(threads, RNGseed=13243223)

  fds <- calculatePSIValues(fds, BPPARAM = bp)

  fds <- filterExpressionAndVariability(
    fds,
    minExpressionInOneSample = minExpressionInOneSample,
    quantileMinExpression = quantileMinExpression,
    quantile = quantile,
    minDeltaPsi = minDeltaPsi,
    filter = filter
  )

  fds <- annotateRangesWithTxDb(
    fds,
    txdb = txdb,
    orgDb = orgdb,
    keytype = "ENTREZID",
    featureName = "SYMBOL"
  )


  fds <- FRASER(fds, BPPARAM = bp)


  return(fds)
}

#' Write FRASER results to file
#'
#' This function extracts FRASER results from a \code{FraserDataSet} object
#' and writes them to a tab-separated file.
#'
#' @param fds A \code{FraserDataSet} object containing fitted FRASER results.
#' @param prefix Character scalar. Prefix for the output filename. The results
#'   are written to \code{<prefix>_fraser_result_full.tsv}.
#' @param test_mode Logical. If \code{TRUE}, indicates test mode (currently
#'   unused). Default: \code{FALSE}.
#'
#' @return Invisibly returns the path to the written file.
#'
#' @details
#' The function calls \code{results(fds, all=TRUE, aggregate=FALSE)} to extract
#' all splice site–level results, and writes them as a tab-separated table
#' without quotes and without row names.
#'
#' @examples
#' \dontrun{
#' fds <- run_fraser(sampleTable, threads = 4, prefix = "fraser_run",
#'                   txdb = txdb, orgdb = orgdb)
#' write_fraser_output(fds, prefix = "fraser_run")
#' }
#'
#' @seealso \code{\link[FRASER]{results}}
#' @export
write_fraser_output <- function(fds, prefix, test_mode=FALSE){

  fraser_full <- results(fds, all=TRUE, aggregate=FALSE)
  write.table(
    fraser_full,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    file = paste0(prefix, "_fraser_result_full.tsv")
  )
}


import_fraser_data <- function(countTable, junctions, splice_sites) {
  sampleTable <- fread(countTable)
  junctions   <- fread(junctions)
  spliceSites <- fread(splice_sites)

  fds <- FraserDataSet(
    colData    = sampleTable,
    junctions  = junctions,
    spliceSites= spliceSites,
    workingDir = "FRASER_output"
  )

  return(fds)
}

#' Run the full FRASER analysis pipeline from command-line arguments
#'
#' This is the main entry point for running the FRASER pipeline. It parses
#' command-line arguments, prepares the sample table and annotation databases,
#' executes the FRASER workflow, and writes results to file.
#'
#' @param args List of command-line arguments (typically from
#'   \code{get_arguments()}). In practice this parameter is ignored, since
#'   \code{get_arguments()} is called inside the function.
#'
#' @return No return value. Side effects include:
#' \itemize{
#'   \item Saving the sample table and annotation objects to
#'     \code{sample_table.RData}.
#'   \item Running FRASER and storing results in the working directory.
#'   \item Writing the FRASER results to a tab-separated file with the given
#'     prefix.
#' }
#'
#' @details
#' The function executes the following steps:
#' \enumerate{
#'   \item Parse arguments with \code{get_arguments()}.
#'   \item Load transcript (\code{TxDb.Hsapiens.UCSC.hg38.knownGene}) and gene
#'     annotation (\code{org.Hs.eg.db}) databases.
#'   \item Build a sample table from input and reference BAM files.
#'   \item Save the sample table and annotation objects to \code{sample_table.RData}.
#'   \item Run the FRASER pipeline via \code{run_fraser()}.
#'   \item Write full FRASER results with \code{write_fraser_output()}.
#' }
#'
#' @examples
#' \dontrun{
#' # Typical usage from command line
#' main()
#' }
#'
#' @seealso \code{\link{get_arguments}}, \code{\link{create_sample_table}},
#'   \code{\link{run_fraser}}, \code{\link{write_fraser_output}}
#' @export
main <- function(args){
  args <- get_arguments()

  taxdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  orgdb <- org.Hs.eg.db

  fds <- import_fraser_data(args$countTable, args$junctionCounts, args$spliceCounts)


  fds <- run_fraser(fds, args$threads, args$prefix, taxdb, orgdb, args$minExpressionInOneSample,
    args$quantileMinExpression, args$quantile, args$minDeltaPsi, args$filter
  )

  write_fraser_output(fds, args$prefix, args$test)
}


main(args)
