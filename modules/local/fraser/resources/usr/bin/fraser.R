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
      metavar="input",
      nargs="+",
      help="Input bam files"
  )
  parser$add_argument(
      "-r", "--refset",
      metavar="refset",
      nargs="+",
      help="Refset bam files"
  )
  parser$add_argument(
      "--prefix",
      metavar="prefix",
      help="Output prefix"
  )
  parser$add_argument(
      "--gtf",
      metavar="gtf",
      help="Gene annotations"
  )
  parser$add_argument(
      "--paired",
      metavar="paired",
      help="Paired-end <true|false>",
      type="logical",
      default=TRUE
  )
  parser$add_argument(
      "-p", "--threads",
      metavar="threads",
      help="Number of parallel threads",
      type="integer",
      default=12
  )
  parser$add_argument(
      "--test",
      metavar="test",
      help="Test run; Just do calculations, don't write outputs"
  )
  parser$add_argument(
      "-v", "--version",
      action="version",
      version=as.character(packageVersion("FRASER"))
  )


  return(parser$parse_args())
}



create_sample_table <- function(sample,ref,paired_end){
  inputs <- c(sample, ref)

  sampleTable <- data.frame(
    sampleID = tools::file_path_sans_ext(basename(inputs)),
    bamFile = inputs,
    stringsAsFactors = FALSE
  )

  sampleTable$pairedEnd <- TRUE



  return(DataFrame(sampleTable))
}


run_fraser <- function(sampleTable, threads, prefix, txdb, orgdb){
  bp <- MulticoreParam(threads, RNGseed=13243223)
  # create FRASER object
  fds <- FraserDataSet(colData=sampleTable, workingDir="./FRASER_output")

  # count reads
  fds <- countRNAData(fds, BPPARAM = bp)

  fds <- calculatePSIValues(fds, BPPARAM = bp)

  fds <- filterExpressionAndVariability(
    fds,
    minExpressionInOneSample = 20,
    quantileMinExpression = 10,
    quantile = 0.95,
    minDeltaPsi = 0.05,
    filter = TRUE
  )



  fds <- annotateRangesWithTxDb(
    fds,
    txdb = txdb,
    orgDb = orgdb,
    keytype = "ENTREZID",
    featureName = "hgncs_symbol"
  )
  fds <- FRASER(fds, BPPARAM = bp)

  return(fds)
}


write_fraser_output <- function(fds, prefix, test_mode=FALSE){

  if (! test_mode ) {
    fraser_full <- results(fds, all=TRUE, aggregate=FALSE)

    write.table(
      fraser_full,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      file = paste0(prefix, "_fraser_result_full.tsv")
    )
  } else if (test_mode){
    # write empty tsv in test mode so nextflow does not error out
    file.create(paste0(prefix, "_fraser_result_full.tsv"))
  }


}


main <- function(args){
  args <- get_arguments()

  taxdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  orgdb <- org.Hs.eg.db
  sampleTable <- create_sample_table(args$input, args$refset, args$paired)
  fds <- run_fraser(sampleTable, args$threads, args$prefix, taxdb, orgdb)
  write_fraser_output(fds, args$prefix)
}


main(args)
