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



args <- commandArgs()

get_arguments <- function(){
  parser <- ArgumentParser(description="Run FRASER")

  parser$add_argument(
    "-i", "--input",
    metavar = "input",
    nargs   = "+",
    help    = "Input bam files"
  )
  parser$add_argument(
    "-r", "--refset",
    metavar = "refset",
    nargs   = "+",
    help    = "Refset bam files"
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



create_sample_table <- function(sample,ref,paired_end){
  inputs <- c(sample, ref)

  Sampletable <- data.frame(
    sampleID = tools::file_path_sans_ext(basename(inputs)),
    bamFile = inputs,
    stringsAsFactors = FALSE
  )

  Sampletable$pairedEnd <- TRUE



  return(DataFrame(Sampletable))
}


run_fraser <- function(sampleTable, threads, prefix, txdb, orgdb, minExpressionInOneSample = 20, quantileMinExpression = 10,
                       quantile = 0.95, minDeltaPsi = 0.05, filter = TRUE){
  bp <- MulticoreParam(threads, RNGseed=13243223)
  # create FRASER object
  fds <- FraserDataSet(colData=sampleTable, workingDir="./FRASER_output")

  # count reads
  fds <- countRNAData(fds, BPPARAM = bp)

  message("Counting done")
  fds <- calculatePSIValues(fds, BPPARAM = bp)

  message("PSI values done")
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

  message("Running fraser")
  fds <- FRASER(fds, BPPARAM = bp)
  message("Fraser itself is done")
  save(fds, file="fraser_data.RData")
  return(fds)
}


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


main <- function(args){
  args <- get_arguments()

  taxdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  orgdb <- org.Hs.eg.db
  sampleTable <- create_sample_table(args$input, args$refset, args$paired)

  save(sampleTable, taxdb, orgdb, file="sample_table.RData")

  fds <- run_fraser(sampleTable, args$threads, args$prefix, taxdb, orgdb, args$minExpressionInOneSample,
    args$quantileMinExpression, args$quantile, args$minDeltaPsi, args$filter
  )

  write_fraser_output(fds, args$prefix, args$test)
}


main(args)
