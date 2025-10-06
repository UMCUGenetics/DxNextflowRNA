#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(FRASER)
  library(argparse)
  library(data.table)
})

 get_arguments <- function(){
  parser <- ArgumentParser(description="Run FRASER")

  parser$add_argument(
    "-r", "--refset",
    metavar = "refset",
    nargs   = "+",
    help    = "Refset bam files"
    )

  parser$add_argument(
    "-c", "--threads",
    metavar = "threads",
    help    = "Number of parallel threads",
    type    = "integer",
    default = 12
  )
  parser$add_argument(
    "-o", "--outdir",
    metavar="outdir",
    default = "./"
  )

  return(parser$parse_args())
}


create_sample_table <- function(ref, paired_end = TRUE){

  Sampletable <- data.frame(
    sampleID = tools::file_path_sans_ext(basename(ref)),
    bamFile = ref,
    stringsAsFactors = FALSE
  )

  Sampletable$pairedEnd <- paired_end

  return(DataFrame(Sampletable))
}

run_fraser <- function(sampleTable, threads){
  bp <- MulticoreParam(threads, RNGseed=13243223)

  fds <- FraserDataSet(colData=sampleTable, workingDir="./FRASER_output")

  fds <- countRNAData(fds, BPPARAM = bp)

  return(fds)
}

export_counts <- function(fds, type, suffix, outdir='./') {

  # retrieve feature metadata
  ranges <- as.data.table(rowRanges(fds, type = type))
  setnames(ranges, old = c("seqnames", "start", "end", "width", "strand"),
                   new = c("seqnames", "start", "end", "width", "strand"))

  # Extract count matrix
  cts <- counts(fds, type = type)

  # Export one file per sample
  for (sample in colnames(cts)) {
    dt <- copy(ranges[, .(seqnames, start, end, width, strand)])
    dt[, sample := as.numeric(cts[, sample])]

    outFile <- file.path(outdir, paste0(sample, suffix))
    fwrite(dt, outFile, sep = "\t", compress = "gzip")

  }
}

main <- function(args){
  args <- get_arguments()

  sampleTable <- create_sample_table(args$refset, paired_end = TRUE)

  fds <- run_fraser(sampleTable,  args$threads)

  # Export junction (j) counts
  export_counts(fds, type="j", suffix="_junction_counts.tsv.gz", args$outdir)

}

main(args)
