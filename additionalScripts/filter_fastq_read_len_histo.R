#!/usr/bin/Rscript

args <- commandArgs(TRUE)
pngFileName <- paste( args[1], ".png", sep="")
data.f <- read.table(args[1], as.is = TRUE, header = TRUE, sep = "\t", check.names=F)




png(filename = pngFileName, width = 800, height = 600, bg = "white")


plot(1, 1, type = "n", xlab = "Genomic location", ylab = "Coverage", axes=F, xlim=c(Xmin, Xmax), ylim=c(0,totalPlotHeight))




dev.off()

