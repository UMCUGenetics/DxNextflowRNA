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
  save.image('./img.RData')

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

## save_output <- function(out_path, out_outrider, ref_samples, prefix, query, padj_thres=0.05, zscore_thres=0, a=TRUE) {
## #    res <- as_tibble(results(out_outrider, padjCutoff=padj_thres, zScoreCutoff=zscore_thres, all=a))
##   ##Reference samples are excluded from final results.
##    #    query_res <- filter(res, !(sampleID %in% ref_samples))
##   res <- as_tibble(results(out_outrider, all=a))
##   query_res <- res

##                                         # Write output table with aberrant expressed targets.
##   write_tsv(query_res, paste0(out_path, prefix, ".outrider_result.tsv"))
## }




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

  ods <- filter_expression(ods, query, prefix, gtf)

                                        # Heatmap of the sample correlation
                                        # it can also annotate the clusters resulting from the dendrogram
  plotCountCorHeatmap(ods, normalized=FALSE)

                                        # Heatmap of the gene/sample expression
  ods <- plotCountGeneSampleHeatmap(ods, normalized=FALSE)

  ods <- estimateSizeFactors(ods)
  message(date(), ": sizeFactors...")
  print(ods$sizeFactor)

  pars_q <- seq(2, min(100, ncol(ods) - 1, nrow(ods) - 1), 2)

  message(date(), ": findEncodingDim...")
  ods <- findEncodingDim(ods, params = pars_q, implementation="autoencoder", BPPARAM=MulticoreParam(nthreads))

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



## ##Only necessary for Erasmus app
## get_ct_emcapp <- function(query, prefix){
##   count_table <- read_delim(query, show_col_types=FALSE, skip=1)
##   ct <- count_table[,c(2,3,4,8,1)]
##   colnames(ct)<- c("chr","start","end","gene","geneID")
##   ct$chr <- str_split_fixed(ct$chr,';',2)[,1]
##   ct$start <- str_split_fixed(ct$start,';',2)[,1]
##   ct$end <- gsub("^.*;", "", ct$end)

##   ##For gene level: create countdata table Identifier chr_start_end_genename
##   ct$emcID <- paste(ct$chr,ct$start,ct$end,ct$gene,sep="_")
##   ##For exon level: create countdata table Identifier: add Ensemble exon_id: chr_start_end_genename-exon_id
##   if(grepl("exon", prefix)){
##     ct$emcID <- paste(ct$emcID,ct$geneID,sep="-")
##   }

##   return(ct)
## }

## ##Only necessary for Erasmus app
## save_res_emcapp <- function(ct, res, prefix, out_path){
##   res <- as_tibble(results(res, all=TRUE))[,c(2,1,3:14)]
##   res_merge <- merge(res, ct, by = "geneID", all.res = TRUE)
##   res_merge$geneENSG <- res_merge$geneID
##   res_merge$geneID <- res_merge$emcID
##   res_merge$sample_name<-NA
##   res_merge$TIN_mean<-NA
##   res_merge$link_bam<-NA

##   fragment <- "genes"
##   treatment <- "untreated"
##   if(grepl("exon",prefix)){ fragment <- "exons"}
##   if(grepl("chx|CHX",prefix)){ treatment <- "CHX"}
##   write_csv(res_merge[,c(1:18,20:23)], paste0(out_path, "umcu_rnaseq_fib_",treatment,"_res_outrider_",fragment,"_counts.tsv"), append=FALSE, col_names = TRUE)
## }

## ##Only necessary for Erasmus app
## save_count_meta_emcapp <- function(ct, all_counts, out_path, prefix){
##   ct <- ct[,c("emcID","geneID")]
##   colnames(all_counts)[1]<-c("geneID")
##   merge_counts <- merge(all_counts, ct, by = "geneID", all.ct = TRUE)
##   merge_counts$geneID <- merge_counts$emcID
##   counts_out <- merge_counts[,c(1:ncol(merge_counts)-1)]
##   colnames(counts_out)[1]<-c("")

##   ##countdata output table EMC app
##   if(grepl("gene", prefix)){
##     write_tsv(counts_out, paste0(out_path, "umcu_rnaseq_genes_counts.tsv"))
##   }else{
##     write_tsv(counts_out, paste0(out_path, "umcu_rnaseq_exons_counts.tsv"))
##   }

##   ##Meta output table EMC app
##   treatment <- "untreated"
##   if(grepl("chx|CHX",prefix)){ treatment <- "chx"}
##   metadata<-data.frame(colnames(counts_out[,2:ncol(counts_out)]),treatment,"fib","umcu_rnaseq",0,"")
##   colnames(metadata)<-c("sample_id","treatment","species","experiment_GS","gender","drop")
##   if(grepl("gene",prefix)){
##     write.csv2(metadata,paste0(out_path, "umcu_rnaseq_metadata_genes.csv"), row.names = FALSE, quote=FALSE)
##   }else{
##     write.csv2(metadata,paste0(out_path, "umcu_rnaseq_metadata_exons.csv"), row.names = FALSE, quote=FALSE)
##   }
## }







main(args$query, args$ref, args$output_path, args$pref, args$gtf, args$nthreads)
