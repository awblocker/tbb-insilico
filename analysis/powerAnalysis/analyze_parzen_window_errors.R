
# Setup -------------------------------------------------------------------


# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

source('powerAnalysis/lib.R')

# Constants

kGroundtruthPath <- 'powerAnalysis/data/simChrom_groundtruth.txt'
kPeakPath <- 'powerAnalysis/parzen_window_peaks.RData'
kReplicates <- 1:10
kPathParzenWindowPower <- 'powerAnalysis/output_parzen_window_power.RData'

# IO ----------------------------------------------------------------------


# Load groundtruth data
groundtruth <- read.csv(kGroundtruthPath)
groundtruth$rep <- groundtruth$rep + 1
groundtruth <- groundtruth[order(groundtruth$rep, groundtruth$pos), ]

# Load Parzen window results
load(kPeakPath)

# Matching ----------------------------------------------------------------


peaks$error <- NA
groundtruth$error <- NA
for (r in kReplicates) {
  # Get indices to for relevant subset of groundtruth, primaries only
  subset.peaks <- which(peaks$chrom == r)
  subset.rep <- which(groundtruth$rep == r & groundtruth$offset == 0)
  
  # Find nearest-neighbor positions
  matched.groundtruth <- MatchPositions(peaks$peaks[subset.peaks], 
                                        groundtruth$pos[subset.rep])
  matched.groundtruth.errors <- peaks$peaks[subset.peaks] - matched.groundtruth
  
  matched.detections <- MatchPositions(groundtruth$pos[subset.rep],
                                       peaks$peaks[subset.peaks])
  matched.detections.errors <- groundtruth$pos[subset.rep] - matched.detections
  
  # Build data frames for detailed analysis
  groundtruth$error[subset.rep] <- matched.detections.errors  
  peaks$error[subset.peaks] <- matched.groundtruth.errors
}

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

# Plotting ----------------------------------------------------------------


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
pdf('powerAnalysis/plots_parzen_window_power.pdf', 10, 6.2, onefile=TRUE)
print(plot.power)
print(plot.error)
dev.off()

# Save results
save(groundtruth, peaks, power.cluster, file=kPathParzenWindowPower)
