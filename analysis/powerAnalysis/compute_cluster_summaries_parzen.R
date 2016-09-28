# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(Rcpp)

source('powerAnalysis/lib.R')
sourceCpp('reproducibility/lib.cpp')

# Constants

kReadsPath <- 'powerAnalysis/data/simChrom_y.txt'
kPeakPath <- 'powerAnalysis/parzen_window_peaks.RData'
kClusterWidth <- 147
kStructureAdj <- 1.
kQuantileSparsity <- 0.9
kOutputPath <- 'powerAnalysis/output_cluster_summaries_parzen.RData'

# IO ----------------------------------------------------------------------

# Load simulated read counts
reads.conn <- file(kReadsPath)
y.line.1 <- readLines(reads.conn, n=1)
y.n.col <- length(strsplit(y.line.1, ',', fixed=TRUE)[[1]])
y.matrix <- matrix(scan(reads.conn, what=integer(0), sep=','),
                   ncol=y.n.col, byrow=TRUE)
close(reads.conn)

# Load Parzen window results
load(kPeakPath)
clusters <- peaks
names(clusters) <- colnames(clusters) <- c('rep', 'center')

# Compute cluster-level summaries
stat.list <- list()
for (replicate in unique(clusters$rep)) {
  subset.clusters.replicate <- which(clusters$rep == replicate)
  
  stat.list[[replicate]] <- ComputeClusterIndices(
    clusters$center[subset.clusters.replicate],
    y.matrix[replicate, ],
    w=as.integer(kClusterWidth %/% 2), adj=kStructureAdj, q=kQuantileSparsity)
  stat.list[[replicate]]$rep <- replicate
}

clusters <- join(clusters, Reduce(rbind, stat.list), type='left')

# Save results
save(clusters, file=kOutputPath)
