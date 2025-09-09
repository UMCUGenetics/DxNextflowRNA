#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(OUTRIDER)
  library(dplyr)
  library(purrr)
  library(magrittr)
  library(readr)
  library(argparse)
  library(BiocParallel)
})

set.seed(3421)


#' Parse command-line arguments for OUTRIDER
#'
#' Defines and parses command-line arguments for running an OUTRIDER analysis.
#' The arguments cover input query/reference files, genome annotation, filtering
#' options, and runtime settings.
#'
#' @return A named list (from \code{argparse::ArgumentParser}) containing the
#'   parsed command-line arguments. Each element corresponds to one of the
#'   options defined below.
#'
#' @section Command-line options:
#' \describe{
#'   \item{\code{-q, --queries}}{One or more query input files.}
#'   \item{\code{-o, --output_path}}{Directory where output will be stored. Default: \code{"./"}.}
#'   \item{\code{-r, --ref}}{One or more reference input files.}
#'   \item{\code{-m, --mode}}{Running mode: \code{"gene"} (default) or \code{"exon"}.}
#'   \item{\code{-g, --gtf}}{Path to the genome annotation in GTF format.}
#'   \item{\code{-t, --threads}}{Number of parallel threads. Default: 12.}
#'   \item{\code{-p, --prefix}}{Output file prefix used for naming results.}
#'   \item{\code{--filter_mode}}{Filtering mode:
#'     \code{"minCounts"} (default) filters only features with zero expression,
#'     \code{"fpkm"} filters additionally lowly expressed features.}
#'   \item{\code{--fpkm_cutoff}}{FPKM cutoff used when \code{--filter_mode fpkm}. Default: 0.5.}
#'   \item{\code{--fpkm_percentile}}{FPKM percentile used when \code{--filter_mode fpkm}. Default: 0.90.}
#'   \item{\code{--mask_samples}}{One or more sample IDs to exclude from model training
#'     (to avoid overfitting dependent samples). Masked samples are still included in
#'     the final results.}
#'   \item{\code{-v, --version}}{Show the OUTRIDER package version and exit.}
#' }
#'
#' @examples
#' \dontrun{
#' args <- get_arguments()
#' print(args$queries)
#' print(args$output_path)
#' }
#'
#' @export
get_arguments <- function() {
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

  parser$add_argument("--fpkm_cutoff", type="double", default=0.5,
                      help="Only used when --filter_mode is fpkm.")

  parser$add_argument("--fpkm_percentile", type="double", default=0.90,
                      help="Only used when --filter_mode is fpkm.")

  parser$add_argument("--mask_samples", metavar="mask_samples", nargs="+",
                      help=paste(c("OUTRIDER expects that each sample is independent.",
                                   "To avoid overfitting in the autoencoder, dependent samples can be masked to exclude ",
                                   "them in the training step. Using this option will keep them in the final output of OUTRIDER")))

  parser$add_argument("-v", "--version", action="version", version=as.character(packageVersion("OUTRIDER")))

  return(parser$parse_args())
}


#' Run OUTRIDER analysis
#'
#' High-level wrapper to perform OUTRIDER aberrant expression/splicing detection
#' using a count matrix and associated sample metadata. The function handles
#' filtering of lowly expressed features, annotation with a GTF, and parallel
#' computation.
#'
#' @param count_matrix \code{matrix} or \code{data.frame}
#'   Raw count matrix with features (e.g., genes, exons, junctions) in rows and
#'   samples in columns. Counts should be unnormalized integer values.
#' @param metadata \code{data.frame}
#'   Sample metadata table containing at least one column with sample IDs
#'   matching the column names of \code{count_matrix}. Can also include covariates
#'   (e.g., batch, condition) used for modeling.
#' @param mode \code{character}
#'   Analysis mode to determine the feature level. Typical options are
#'   \code{"gene"} for gene-level counts, \code{"exon"} for exon usage, or
#'   \code{"junction"} for splice-junction counts.
#' @param gtf \code{character}
#'   File path to a gene annotation in GTF format. Used to map features to genes
#'   and provide genomic context.
#' @param filter_mode \code{character}
#'   Filtering strategy for removing lowly expressed features. Common options are
#'   \code{"fpkm"} or \code{"percentile"}, depending on whether absolute or
#'   relative expression thresholds are applied.
#' @param fpkm_cutoff \code{numeric}
#'   Minimum FPKM threshold used when \code{filter_mode = "fpkm"}. Features with
#'   FPKM values below this cutoff across all samples are discarded.
#' @param fpkm_percentile \code{numeric}
#'   Percentile threshold (0–100) used when \code{filter_mode = "percentile"}.
#'   Features below this expression percentile are removed.
#' @param threads \code{integer}
#'   Number of parallel threads to use during computation. Increasing this value
#'   speeds up runtime on multicore systems.
#'
#' @return An \code{OutriderDataSet} object containing filtered counts, fitted
#'   model parameters, quality control metrics, and significance results for
#'   aberrant expression/splicing.
#'
#' @examples
#' \dontrun{
#' result <- run_outrider(
#'     count_matrix     = counts,
#'     metadata         = sample_info,
#'     mode             = "gene",
#'     gtf              = "annotations.gtf",
#'     filter_mode      = "fpkm",
#'     fpkm_cutoff      = 1,
#'     fpkm_percentile  = 95,
#'     threads          = 8
#' )
#' }
#'
run_outrider <- function(count_matrix, metadata, mode, gtf, filter_mode, fpkm_cutoff, fpkm_percentile, threads){
  bp <- MulticoreParam(threads, RNGseed = 13243223)
  ods <- OutriderDataSet(countData = count_matrix)

  if(mode == "gene"){
      if( filter_mode == "minCounts"){
          ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE)
      } else {
        ods <- filterExpression(
          ods,
          fpkmCutoff = fpkm_cutoff,
          percentile = fpkm_percentile,
          minCounts = FALSE,
          filterGenes = TRUE,
          gtfFile = gtf)
      }

  } else if(mode == "exon") {
    # Add exon length to the outrider dataset to ensure FPKM is correctly estimated
    mcols(ods)$basepairs <- metadata$Length
    ods <- filterExpression(
      ods,
      fpkmCutoff = fpkm_cutoff,
      percentile  = fpkm_percentile,
      filterGenes = TRUE
    )
  }

  # Runs the outrider wrapper in parallel
  ods <- OUTRIDER(
    ods,
    BPPARAM=bp
  )

  return(ods)
}

#' Read and split count table
#'
#' Reads a tab-delimited count table and splits it into metadata columns and a
#' count column. This function assumes that the last column always contains raw
#' counts, while all preceding columns represent feature metadata.
#'
#' @param file \code{character}
#'   Path to a tab-delimited text file with a header. The file must contain one
#'   or more metadata columns followed by a final column with integer counts.
#'
#' @return A \code{list} with two elements:
#'   \item{metadata}{\code{data.frame} with all columns except the last.}
#'   \item{counts}{\code{numeric} vector with the values from the last column.}
#'
#' @examples
#' \dontrun{
#' result <- read_and_split_counts("counts_with_metadata.tsv")
#' head(result$metadata)
#' head(result$counts)
#' }
#'
read_and_split_counts <- function(file) {
    df <- read.table(file, sep="\t", header=T)
    count_col <- df[[ncol(df)]]
    meta_cols <- df[, -ncol(df)]
    list(metadata = meta_cols, counts = count_col)
}

#' Main entry point for OUTRIDER analysis
#'
#' This function coordinates the full OUTRIDER workflow. It:
#' \enumerate{
#'   \item Parses command-line arguments via \code{\link{get_arguments}}.
#'   \item Reads query and reference count tables using \code{\link{read_and_split_counts}}.
#'   \item Validates that metadata columns are consistent across all input files.
#'   \item Combines individual count vectors into a joint count matrix.
#'   \item Runs the OUTRIDER model via \code{\link{run_outrider}}.
#'   \item Collects both full and significant results.
#'   \item Writes results to tab-delimited files with the specified prefix.
#' }
#'
#' @return Invisibly returns \code{NULL}. Side effects include reading input
#'   files, running the OUTRIDER analysis, and writing two result tables to disk:
#'   \itemize{
#'     \item \code{<prefix>_outrider_result_full.tsv}: full results without FDR filtering
#'     \item \code{<prefix>_outrider_result_signif.tsv}: significant results (after FDR filtering)
#'   }
#'
#' @seealso \code{\link{get_arguments}}, \code{\link{read_and_split_counts}},
#'   \code{\link{run_outrider}}
#'
#' @examples
#' \dontrun{
#' main()
#' }
#'
main <- function(){

  args <- get_arguments()


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

  ods <- run_outrider(count_matrix, metadata, args$mode, args$gtf, args$filter_mode, args$fpkm_cutoff, args$fpkm_percentile, 
                      args$threads)
  res_full <- results(ods, all = TRUE)
  res_signif <- results(ods, all = FALSE)

  # Output both the unfiltered results, and the filtered results
  write.table(
    res_full,
    quote = F, sep = "\t", row.names = FALSE,
    file = paste0(args$prefix,"_outrider_result_full.tsv")
  )
  write.table(
    res_signif,
    quote = F, sep = "\t", row.names = FALSE,
    file = paste0(args$prefix,"_outrider_result_signif.tsv")
  )

}



main()
