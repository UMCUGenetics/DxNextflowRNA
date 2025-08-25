#!/usr/bin/env Rscript

library(FRASER)
library(biomaRt)
library(utils)
library(argparse)

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

    message(date(), ": sampleTable")
    print(sampleTable)

    return(DataFrame(sampleTable))
}


run_fraser <- function(sampleTable, threads, prefix){
    # create FRASER object
    fds <- FraserDataSet(colData=sampleTable, workingDir="./FRASER_output")

    print(fds)

    ## register(MulticoreParam(workers=min(threads, multicoreWorkers())))

    # count reads
    fds <- countRNAData(fds, BPPARAM = MulticoreParam(threads))



    message(date(), ": RNA counts calulated")
    fds <- calculatePSIValues(fds)

    message(date(), ": PSI values calulated")
    ## FRASER vignette/gagneurlab settings:
    ## Currently, we suggest keeping only junctions which support the following:
    ## minExpressionInOneSample - 20 - At least one sample has 20 (or more) reads
    ## quantileMinExpression - 10 - The minimum total read count (N) an intron needs to have at the specified quantile across samples to pass the filter. See quantileForFiltering.
    ## quantile - 0.95 -  Defines at which percentile the quantileMinExpression filter is applied.  # nolint
    ## A value of 0.95 means that at least 5% of the samples need to have a total read count N >= quantileMinExpression to pass the filter.
    ## minDeltaPsi - 0.05 -  The minimal variation (in delta psi) required for an intron to pass the filter. gagneurlab suggests 0.05

  fds <- filterExpressionAndVariability(fds,  minExpressionInOneSample = 20, quantileMinExpression=10, quantile = 0.95, minDeltaPsi=0.05, filter=FALSE)
  message(date(), ": filter done")
    fds <- saveFraserDataSet(fds, name="fds-filter")


    pdf(paste0(prefix,"_heatmap.pdf"))
    print(plotFilterExpression(fds, bins=100))
    dev.off()

    message(date(), ": dimension before filtering...")
    print(dim(assays(fds)$rawCountsJ))
    print(dim(assays(fds)$rawCountsSS))

    fds <- fds[mcols(fds, type="j")[,"passed"]]
    message(date(), ": dimension after filtering...")
    print(dim(assays(fds)$rawCountsJ))
    print(dim(assays(fds)$rawCountsSS))

    #Finds the optimal encoding dimension by injecting artificial splicing outlier ratios while maximizing the precision-recall curve.
    # Get range for latent space dimension, max tested dimensions = 6
    mp <- 6
    a <- 2
    b <- min(ncol(fds), nrow(fds)) / mp   # N/mp

    maxSteps <- 12
    if(mp < 6){
        maxSteps <- 15
    }

    Nsteps <- min(maxSteps, b)
    pars_q <- unique(round(exp(seq(log(a),log(b),length.out = Nsteps))))

    message(date(), ": getEncDimRange", pars_q)

    fds <- saveFraserDataSet(fds, name="fds-normfalse")
    # Heatmap of the sample correlation
    # plotCountCorHeatmap(fds, type="jaccard", logit=FALSE, normalized=FALSE)
    # Heatmap of the intron/sample expression
    # plotCountCorHeatmap(fds, type="jaccard", logit=FALSE, normalized=FALSE, plotType="junctionSample", topJ=100, minDeltaPsi=0.05)

    #set.seed(as.integer(random))
    # hyperparameter opimization can be used to estimate the dimension q of the latent space of the data.
    # It works by artificially injecting outliers into the data and then comparing the AUC of recalling these outliers for different values of q.
    fds <- optimHyperParams(fds, type="jaccard", implementation="PCA",
    q_param=pars_q, minDeltaPsi=0.05, plot=FALSE)
    # retrieve the estimated optimal dimension of the latent space
    bestq = bestQ(fds, type="jaccard")
    message(date(), "...find best q: " ,bestq)

    # Add verbosity to the FRASER object
    verbose(fds) <- 3
    fds <- fit(fds, q=bestq, type="jaccard", implementation="PCA", iterations=15)

    # Heatmap of the sample correlation
    #  plotCountCorHeatmap(fds, type="jaccard", logit=FALSE, normalized=TRUE)
    # Heatmap of the intron/sample expression
    #  plotCountCorHeatmap(fds, type="jaccard", logit=FALSE, normalized=TRUE, plotType="junctionSample", topJ=100, minDeltaPsi=0.05)

    fds <- saveFraserDataSet(fds, name="fds-normtrue")

    #all in one getEncDimRange, optimHyperParams, bestQ, fit
    #fds_fraser <- FRASER(fds_filtered, q=c(jaccard=2))

    fds <- annotateRanges(fds, GRCh = 38)

    # Pvalues
    fds <- calculatePvalues(fds, type="jaccard")
    # Adjust Pvalues
    fds <- calculatePadjValues(fds, type="jaccard")

    fds <- saveFraserDataSet(fds, name="fds-final")

    return(fds)
}


write_fraser_output <- function(fds, prefix){
    # retrieve results with default and recommended cutoffs (padj <= 0.05 (0.1 gagneurlab) and # |deltaPsi| >= 0.3)
    fraser_res <- results(fds, all=TRUE)
    message(date(), ": result...")
    print(fraser_res)

    write.table(as.data.frame(unname(fraser_res)), paste0("fraser_result_",prefix,".tsv"), sep="\t", quote=FALSE)
}


main <- function(args){
    args <- get_arguments()

    sampleTable <- create_sample_table(args$input, args$refset, args$paired)
    fds <- run_fraser(sampleTable, args$threads, args$prefix)
    write_fraser_output(fds, args$prefix)


}


main(args)
