library(plyr)

# Set clusters
kClusterPattern <- "results/mcmc_clusters_%s_chrom%02d.txt"
kDetectionPattern <- "results/mcmc_detections_%s_chrom%02d_pm%d.txt"
kExperiments <- c('H_1-1', 'H_1-2', 'H_2-1', 'H_2-2',
                  'H_1-combined', 'H_2-combined',
                  sprintf('titration-%d', 1:4))
kChromosomes <- 1:16
kPm <- 0:3


# Script

# Load clusters
data.list <- expand.grid(kExperiments, kChromosomes)
names(data.list) <- c('experiment', 'chrom')

clusters <- mdply(data.list, function(experiment, chrom, pattern) {
  read.table(sprintf(pattern, experiment, chrom), header=TRUE)
}, pattern=kClusterPattern)
                                                                
# Load detections
data.list <- expand.grid(kExperiments, kChromosomes, kPm)
names(data.list) <- c('experiment', 'chrom', 'pm')

detections <- mdply(data.list, function(experiment, chrom, pm, pattern) {
  read.table(sprintf(pattern, experiment, chrom, pm), header=TRUE)
}, pattern=kDetectionPattern)

# Save results
save(clusters, detections, file='out/results.RData')
