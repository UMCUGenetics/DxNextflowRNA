#! /usr/bin/env Rscript
#Rscript MAE.r input_file(tsv file from GATK)  output_file(tsv file)

#import required libraries
library(tMAE)
library(ggplot2)

args<-commandArgs(TRUE)

#import data
countdata <- read.table(args[1], header=TRUE, sep="\t")

# run MAE analysis

result_MAE <- DESeq4MAE(countdata, minCoverage = 10)
result_MAE[, signif := padj < 0.05]
result_MAE[signif == TRUE, .N]
result_MAE[, signif_ALT := signif == TRUE & altRatio >= 0.8]
result_MAE[signif_ALT == TRUE, .N]

#save the pdf plot
pdf("data distribution plot.pdf")
ggplot(countdata, aes(refCount+1, altCount+1)) +
        geom_point() +
        geom_abline(slope=1, intercept=0) +
        scale_y_log10() +
        scale_x_log10() +
        theme_bw()
dev.off()

#export the MAE in tsv file
write.table(result_MAE, file = args[2], sep = "\t")
