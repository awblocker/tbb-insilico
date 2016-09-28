# Load libraries
library(plyr)
library(stringr)
library(reshape2)
library(Rcpp)
library(RcppCNPy)

source('powerAnalysis/lib.R')
sourceCpp('reproducibility/lib.cpp')

# Constants
kRegionPath <- 'powerAnalysis/data/simChrom_regions.txt'  # Path for regions
kCoefPath <- 'powerAnalysis/data/simChrom_b.npy'  # Path for coefficients
kGeneLength <- 3501
kPromoterLength <- 1000
kGroundtruthPath <- 'powerAnalysis/data/simChrom_groundtruth.txt'
kReplicates <- 1:10
kClusterWidth <- 147
kStructureAdj <- 0.
kQuantileSparsity <- 0.9
kOutputPath <- 'powerAnalysis/output_cluster_summaries_groundtruth.RData'


# Script


# Load coefficients
b.matrix <- npyLoad(kCoefPath)

# Load groundtruth data
groundtruth <- read.csv(kGroundtruthPath)
groundtruth$rep <- groundtruth$rep + 1
groundtruth <- groundtruth[order(groundtruth$rep, groundtruth$pos), ]

# Subset to primaries
primaries <- groundtruth[groundtruth$type==0, ]
clusters <- primaries[, c('rep', 'pos', 'coverage', 'eff.magnitude',
                          'magnitude', 'offset')]
names(clusters)[2] <- colnames(clusters)[2] <- 'center'
clusters <- clusters[order(clusters$rep, clusters$center), ]

# Compute cluster-level summaries
stat.list <- list()
for (replicate in unique(clusters$rep)) {
  subset.clusters.replicate <- which(clusters$rep == replicate)
  
  stat.list[[replicate]] <- ComputeClusterIndices(
    clusters$center[subset.clusters.replicate],
    b[replicate, ],
    w=as.integer(kClusterWidth %/% 2), adj=kStructureAdj, q=kQuantileSparsity)
  stat.list[[replicate]]$rep <- replicate
}

clusters <- join(clusters, Reduce(rbind, stat.list), type='left')

save(clusters, file=kOutputPath)
