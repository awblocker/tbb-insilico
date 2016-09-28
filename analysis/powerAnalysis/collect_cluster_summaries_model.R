# Load libraries
library(plyr)
library(stringr)
library(reshape2)

source('powerAnalysis/lib.R')

# Constants

kClustersPattern <- 'results/mcmc_clusters_%s_chrom%02d.txt'
kReplicates <- 1:10
kClustersFmt <- list(center=integer(0),
                     cluster_length=integer(0),
                     occupancy=double(0),
                     occupancy_se=double(0),
                     localization=double(0),
                     localization_se=double(0),
                     structure=double(0),
                     structure_se=double(0),
                     sparsity=double(0),
                     sparsity_se=double(0)
                     )
kOutputPath <- 'powerAnalysis/output_cluster_summaries_model.RData'


# Script

# Load data from all replicates
clusters <- ldply(kReplicates, function(rep)
  data.frame(rep=rep, scan(sprintf(kClustersPattern, 'powerAnalysis', rep),
                           what=kClustersFmt, skip=1),
             rep=rep))

save(clusters, file=kOutputPath)
