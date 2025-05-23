#!/usr/bin/env Rscript

library("argparse")

library("OUTRIDER")
library("tools")
library("tibble")
library("readr")
library("dplyr")
library("stringr")


parser <- ArgumentParser(description = "Process some integers")
parser$add_argument(
         "query",
         metavar = "query_input_files",
         nargs = "+",
         help = "Files or directories containing the query input files."
       )
parser$add_argument(
         "-o",
         "--output_path",
         metavar="output_path",
         nargs="+",
         help="Path where output of OUTRIDER will be stored.",
         default="./"
       )
parser$add_argument(
         "-r",
         "--ref",
         metavar="reference_input_files",
         nargs="+",
         help="Files or directories containing the reference input files."
       )
parser$add_argument(
         "-p",
         "--pref",
         metavar="prefix",
         nargs="+",
         help="Prefix of file.",
         default="gene"
       )
parser$add_argument(
         "-g",
         "--gtf",
         metavar="genome_gtf",
         nargs="+",
         help="Genome gtf file"
       )
parser$add_argument(
         "-t",
         "--threads",
         metavar="nthreads",
         default=4,
         help="Number of parallel threads",
         type="integer"
       )

args <- parser$parse_args()




main <- function(query, ref, output_path, prefix, gtf, nthreads){


  query_data <- get_input(query)
  ref_data   <- get_input(ref)

  all_counts <- merge_count_tables(query_data$count_tables, ref_data$count_tables)
  all_counts.df <- as.data.frame(all_counts)

  rownames(all_counts.df) <- all_counts.df$rownames
  all_counts.df$rownames <- NULL
  all_counts.df$gene_name.x <- NULL
  all_counts.df$gene_name.y <- NULL


  ods <- OutriderDataSet(countData = all_counts.df)
  ## ods <- filterExpression(ods, gtf)
  ods <- filter_expression(ods, query, prefix, gtf)

  ods <- run_outrider(ods, query, prefix, gtf, nthreads)

  res <- results(ods)

  write_tsv(query_res, paste0(out_path, prefix, ".outrider_result.tsv"))
}




read_input_files <- function(input){
  sampleIDs <- c()
  count_tables <- lapply(input, function(f) {
    input_ext <- tools::file_ext(f)
    if(input_ext == "Rds") {
      rse <- readRDS(f)  # RangedSummarizedExperiment
      sampleIDs <- append(sampleIDs, colnames(rse))
      return(as_tibble(rownames_to_column(as.data.frame(assays(rse)$counts), var="rownames")))
    } else if (input_ext %in% c("txt", "tsv")) {
      count_table <- read_delim(f, show_col_types=FALSE, skip=1)
                                        # col 1: samples IDs, col 9: counts
      ct <- count_table[,c(1,9)]
      names(ct)[1] <- "rownames"
      sampleIDs <<- append(sampleIDs, colnames(ct)[2])
      return(ct)
    } else {
      stop("Input file extension is not supported.")
    }
  })
  return(list("sampleIDs"=sampleIDs, "count_tables"=count_tables))
}


get_input <- function(input){
                                        # If directories are provided
  if(all(sapply(input, function(x) dir.exists(x)))) {
    retrieved_files <- sapply(input, function(d){
      list.files(path = d, pattern = "Rds|txt|tsv", full.names = TRUE, recursive = FALSE)
    })
  } else if(all(sapply(input, function(x) file.exists(x)))) {  # If files are provided.
    retrieved_files <- input
  } else {  # If both dir and files are provided, or different tyes (character, int etc)
    stop("Input is neither dir or file.")
  }
  return(read_input_files(retrieved_files))
}

merge_count_tables <- function(lst_query, lst_ref){
                                        # merge count tables together.
  lst_count_tables <- c(lst_query, lst_ref)
  all_counts <- lst_count_tables %>%
    Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1, dtf2,by="rownames"), .)
  return(all_counts)
}

run_outrider <- function(ods, query, prefix, gtf, nthreads) {

  plotheat = "counts_heatplots.pdf"
  pdf(plotheat,onefile = TRUE)

                                        # Heatmap of the sample correlation
                                        # it can also annotate the clusters resulting from the dendrogram
  plotCountCorHeatmap(ods, normalized=FALSE)

                                        # Heatmap of the gene/sample expression
  ods <- plotCountGeneSampleHeatmap(ods, normalized=FALSE)

  ods <- estimateSizeFactors(ods)
  message(date(), ": sizeFactors...")
  print(ods$sizeFactor)

  ## pars_q <- seq(2, min(100, ncol(ods) - 1, nrow(ods) - 1), 2)

  ## message(date(), paste(c(": ", pars_q), collapse=", "))
  message(date(), ": findEncodingDim...")
  ods <- findEncodingDim(ods,
                         ## params=pars_q,
                         implementation="autoencoder",
                         BPPARAM=MulticoreParam(nthreads)
                         )

                                        #find best q (dimension)
  opt_q <- getBestQ(ods)
  message(date(), ": getbestQ...")
  print(opt_q)

  ods <- controlForConfounders(ods, q = opt_q, BPPARAM = MulticoreParam(nthreads))

                                        # After controlling for confounders the heatmap should be plotted again.
                                        # If it worked, no batches should be present and the correlations between samples should be reduced and close to zero. [1]*
  ods <- plotCountCorHeatmap(ods, normalized=TRUE)
  ods <- plotCountGeneSampleHeatmap(ods, normalized=TRUE)
  dev.off()

                                        #    if(grepl("^(peer|pca)$", implementation)){
                                        #         message(date(), ": Fitting the data ...") # NOT with autoencoder
                                        #         ods <- fit(ods, BPPARAM=MulticoreParam(8))
                                        #    }

  ods <- computePvalues(ods, alternative="two.sided", method="BY", BPPARAM = MulticoreParam(nthreads))

  ods <- computeZscores(ods)

                                        #    run full OUTRIDER pipeline (control, fit model, calculate P-values)
                                        #    out <- OUTRIDER(ods, BPPARAM=MulticoreParam(8))
  return(ods)
}



filter_expression <- function(ods, query, prefix, gtf){
  if(grepl("gene", prefix)){
    ##FOR GENE LEVEL
    ods <- filterExpression(ods, fpkmCutoff = 1, minCounts = FALSE, filterGenes = FALSE, gtfFile=gtf)
  }
  else{
    ##FOR EXON LEVEL
    ct <- read_delim(query, show_col_types=FALSE, skip=1)
    print("identical exon ids:")
    print(identical(ct$Geneid,rownames(assays(ods)$counts)))
    readD <- apply(assays(ods)$counts, 2, function(x) x / sum(x) * 10^6)
    countsFPKM <- readD / ct$Length * 10^3
    perc95e <- apply(countsFPKM, 1, function(x) quantile(x,probs=0.95))
    mcols(ods)$passedFilter <- perc95e>1
    mcols(ods)$basepairs <- ct$Length

    ##FOR EXON RATIO
    ##REMOVE ##FOR EXON LEVEL code above
    ##UNCOMMENT BELOW STEPS
                                        #       print("EXON RATIO")
                                        #       ods <- filterExpression(ods, minCounts = TRUE, filterGenes = FALSE)
  }

                                        # display the FPKM distribution of counts.
                                        #ods <- plotFPKM(ods)
  message(date(), ": dim before filtering...")
  print(dim(assays(ods)$counts))

  ##SUBSETTING
  ods <- ods[mcols(ods)$passedFilter,]
  message(date(), ": dim after filtering...")
  print(dim(assays(ods)$counts))
  return(ods)
}



main(args$query, args$ref, args$output_path, args$pref, args$gtf, args$threads)
