# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

source('powerAnalysis/lib.R')

# Constants

kGroundtruthPath <- 'powerAnalysis/data/simChrom_groundtruth.txt'
kClusterPattern <- 'results/mcmc_clusters_powerAnalysis_chrom%02d.txt'
kReplicates <- 1:10
kPathClusterPower <- 'powerAnalysis/output_cluster_power.RData'


# Load groundtruth data
groundtruth <- read.csv(kGroundtruthPath)
groundtruth$rep <- groundtruth$rep + 1
groundtruth <- groundtruth[order(groundtruth$rep, groundtruth$pos), ]

clusters.list <- list()
groundtruth$error <- NA
for (r in kReplicates) {
  # Get indices to for relevant subset of groundtruth
  subset.rep <- which(groundtruth$rep == r & groundtruth$offset == 0)
  
  # Load detections for each replicate
  clusters.list[[r]] <- read.table(sprintf(kClusterPattern, r), header=TRUE)
  clusters.list[[r]]$rep <- r
  
  # Find nearest-neighbor positions
  matched.groundtruth <- MatchPositions(clusters.list[[r]]$center, 
                                        groundtruth$pos[subset.rep])
  matched.groundtruth.errors <- clusters.list[[r]]$center - matched.groundtruth
  
  matched.detections <- MatchPositions(groundtruth$pos[subset.rep],
                                       clusters.list[[r]]$center)
  matched.detections.errors <- groundtruth$pos[subset.rep] - matched.detections
  
  # Build data frames for detailed analysis
  groundtruth$error[subset.rep] <- matched.detections.errors  
  clusters.list[[r]]$error <- matched.groundtruth.errors
}

# Combined detections information
clusters <- Reduce(rbind, clusters.list)
rm(clusters.list)
gc()

# Label primaries with design offsets
groundtruth$offset.primary <- 0
is.na(groundtruth$offset.primary) <- groundtruth$type > 0
ind <- which(groundtruth$type == 0 & groundtruth$eff.magnitude < 1)
groundtruth$offset.primary[ind] <- groundtruth$offset[ind + 1]

# Compute error indicators
groundtruth$detected.primary <- (abs(groundtruth$error) < 5)

power.cluster <- ddply(groundtruth[groundtruth$type==0, ],
                       .(magnitude, coverage, eff.magnitude, offset.primary),
                       summarise,
                       mean.err=mean(abs(error)),
                       sd.err=sd(abs(error)),
                       se.err=sd(abs(error))/sqrt(length(error)),
                       median.err=median(abs(error)),
                       mean.sq.err=mean(error^2),
                       p.detected.primary=mean(detected.primary),
                       n=length(detected.primary))
power.cluster$coverage.bin <- as.integer(factor(power.cluster$coverage))
power.cluster$eff.magnitude <- round(power.cluster$eff.magnitude, 3)

# Build basic summary plots
plot.power <- qplot(x=coverage, y=p.detected.primary,
                    data=power.cluster,
                    geom=c('line', 'point'), ylim=c(0,1))
plot.power <- plot.power + facet_grid(
  offset.primary ~ eff.magnitude, labeller=function(var, value) round(value,3))
plot.power <- plot.power + labs(list(
    x='Coverage',
    y='Detection probability',
    title='Cluster center detection probability by offset and effective magnitude'))

plot.error <- qplot(x=coverage, y=mean.err,
                    ymin=mean.err - 1.96 * se.err,
                    ymax=mean.err + 1.96 * se.err,
                    data=power.cluster,
                    geom=c('line', 'errorbar'))
plot.error <- plot.error + facet_grid(
  offset.primary ~ eff.magnitude, labeller=function(var, value) round(value,3))
plot.error <- plot.error + labs(list(
  x='Coverage',
  y='Mean error',
  title='Cluster center errors by offset and effective magnitude'))


# Save plots
pdf('powerAnalysis/plots_cluster_power.pdf', 10, 6.2, onefile=TRUE)
print(plot.power)
print(plot.error)
dev.off()

# Save results
save(groundtruth, clusters, power.cluster, file=kPathClusterPower)
