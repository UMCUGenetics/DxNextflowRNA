#! /usr/bin/env Rscript
#Rscript MAE.r input_file output_file

#import required libraries
library(tMAE)
library(ggplot2)

args<-commandArgs(TRUE)



#import data
mydata <- read.table(args[1], header=TRUE, sep="\t")


# run MAE analysis

		final <- resMAE <- DESeq4MAE(mydata, minCoverage = 10)
	resMAE[, signif := padj < 0.05]
resMAE[signif == TRUE, .N]
	resMAE[, signif_ALT := signif == TRUE & altRatio >= 0.8]
		resMAE[signif_ALT == TRUE, .N]

#plot the data
ggplot(mydata, aes(refCount+1, altCount+1)) + geom_point() +
    geom_abline(slope=1, intercept=0) + 
    scale_y_log10() + scale_x_log10() + theme_bw()

#export the MAE in csv file
write.table(final, file = args[2], sep = "\t")
