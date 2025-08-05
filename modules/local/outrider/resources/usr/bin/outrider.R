#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(OUTRIDER))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(BiocParallel))



parser <- ArgumentParser(description = "Run OUTRIDER")
parser$add_argument("-q", "--queries", metavar="query_files", nargs="+",
                    help="Files containing the query input files.")

parser$add_argument("-o", "--output_path", metavar="output_path", default="./",
                    help="Path where output will be stored.")

parser$add_argument("-r", "--ref", metavar="reference_files", nargs="+",
                    help="Files containing the reference input files.")

parser$add_argument("-m", "--mode", metavar="mode", default="gene",
                    help="Running mode, gene or exon")

parser$add_argument("-g", "--gtf", metavar="genome_gtf",
                    help="Genome gtf file")

parser$add_argument( "-t", "--threads", metavar="threads", default=12, type="integer",
                    help="Number of parallel threads")

parser$add_argument("-p", "--prefix", metavar="prefix",
                    help="Number of parallel threads")

parser$add_argument("--filter_mode", choices=c("minCounts","fpkm"), default="minCounts",
                    help=paste(c("Gene/Exon filter mode. minCounts filters only genes/exons with zero expression,",
                                 "whereas fpkg filters also lowly epxressed genes/exons")))

parser$add_argument("--fpkm_cutoff", type="double", default=1.0,
                    help="Only used when --filter_mode is fpkm.")

parser$add_argument("--fpkm_percentile", type="double", default=0.95,
                    help="Only used when --filter_mode is fpkm.")


parser$add_argument("-v", "--version", action="version", version=as.character(packageVersion("OUTRIDER")))

args <- parser$parse_args()
set.seed(1)
bp <- MulticoreParam(args$threads, RNGseed=13243223)

run_outrider <- function(count_matrix, metadata, args){

  ods <- OutriderDataSet(countData = count_matrix)

  if(args$mode == "gene"){
      if( args$filter_mode == "minCounts"){
          ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
      } else {
        ods <- filterExpression(
          ods,
          fpkmCutoff = args$fpkm_cutoff,
          percentile = args$fpkm_percentile,
          minCounts = FALSE,
          filterGenes = TRUE,
          gtfFile = args$gtf)
      }

  } else if(args$mode == "exon") {
    mcols(ods)$basepairs <- metadata$Length # add exon length to the outrider dataset to ensure FPKM is correctly estimated
    ods <- filterExpression(
      ods,
      fpkmCutoff = args$fpkm_cutoff,
      percentile  = args$fpkm_percentile,
      filterGenes = TRUE
    )
  }

  # Runs the outrider pipeline in parallel
  ods <- OUTRIDER(
    ods,
    BPPARAM=bp
  )

  return(ods)
}

read_and_split_counts <- function(file) {
    df <- read.table(file, sep="\t", header=T)
    count_col <- df[[ncol(df)]]
    meta_cols <- df[, -ncol(df)]
    list(metadata = meta_cols, counts = count_col)
}

main <- function(args){

  count_files <- c(args$queries, args$ref)

  all_data <- map(count_files, read_and_split_counts)
  all_metadata <- map(all_data, "metadata")

  # Check if metadata is the same between samples
  # If not the genes/exons are not the same and the counts will not match
  if (!all(map_lgl(all_metadata, ~ all.equal(.x, all_metadata[[1]]) == TRUE))) {
    stop("Metadata columns differ across files")
  }

  metadata <- all_metadata[[1]]

  # Combine sample specific count data.frames into a joint count matrix
  count_matrix <- map_dfc(all_data, "counts") %>%
    set_names(basename(count_files)) %>%
    as.matrix()

  rownames(count_matrix) <- metadata$Geneid

  ods <- run_outrider(count_matrix, metadata, args)
  res_full <- results(ods, all=TRUE)
  res_signif <- results(ods, all=FALSE)

  # Output both the unfiltered results, as the filtered results
  write.table(
    res_full,
    quote=F, sep="\t", row.names=F,
    file=paste0(args$prefix,"_outrider_result_full.tsv")
  )
  write.table(
    res_signif,
    quote=F, sep="\t", row.names=F,
    file=paste0(args$prefix,"_outrider_result_signif.tsv")
  )

}



main(args)
