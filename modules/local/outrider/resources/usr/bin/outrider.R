#!/usr/bin/env Rscript

library(OUTRIDER)
library(dplyr)
library(purrr)
library(magrittr)
library(readr)
library(argparse)

parser <- ArgumentParser(description = "Run OUTRIDER")
parser$add_argument(
    "-q", "--queries",
    metavar="query_files",
    nargs="+",
    help="Files containing the query input files."
)
parser$add_argument(
    "-o", "--output_path",
    metavar="output_path",
    help="Path where output will be stored.",
    default="./"
)
parser$add_argument(
    "-r", "--ref",
    metavar="reference_files",
    nargs="+",
    help="Files containing the reference input files."
)
parser$add_argument(
    "-m", "--mode",
    metavar="mode",
    help="Running mode, gene or exon",
    default="gene"
)
parser$add_argument(
    "-g", "--gtf",
    metavar="genome_gtf",
    help="Genome gtf file"
)
parser$add_argument(
    "-t", "--threads",
    metavar="threads",
    default=12,
    help="Number of parallel threads",
    type="integer"
    )
parser$add_argument(
    "-p", "--prefix",
    metavar="prefix",
    help="Number of parallel threads"
)
parser$add_argument(
  "-v", "--version",
  action="version",
  version=as.character(packageVersion("OUTRIDER"))
)

args <- parser$parse_args()


run_outrider <- function(count_matrix, metadata, mode, gtf, threads){

  ods <- OutriderDataSet(countData = count_matrix)

  if(mode == "gene"){
    ods <- filterExpression(
      ods,
      fpkmCutoff=1,
      minCounts=TRUE,
      filterGenes=TRUE,
      gtfFile=gtf
    )
  } else if(mode == "exon") {
    mcols(ods)$basepairs <- metadata$Length # add exon length to the outrider dataset to ensure FPKM is correctly estimated
    ods <- filterExpression(
      ods,
      fpkmCutoff = 0.5,
      percentile  = 0.90,
      filterGenes = TRUE
    )
  }

  ods <- OUTRIDER(
    ods,
    BPPARAM=MulticoreParam(threads)
  )

  return(ods)
}

make_query_plots <- function(query_files, ods, prefix, outdir = "./"){
  volcano_dir <- paste0(outdir, prefix, "_volcano_plots/")
  dir.create(volcano_dir, showWarnings = FALSE)

  for (query in query_files){
    file_name <- paste0(volcano_dir, query, ".png")
    message(paste("Plotting volcano for: ", file_name, "In:", volcano_dir ))
    pdf(file_name)
    print(plotVolcano(ods, sampleID=query, basePlot=T))
    dev.off()
  }
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

  ods <- run_outrider(count_matrix, metadata, mode=args$mode, gtf=args$gtf, threads=args$threads)
  res_full <- results(ods, all=TRUE)
  res_signif <- results(ods, all=FALSE)

  make_query_plots(args$queries, ods, prefix=args$prefix)

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
