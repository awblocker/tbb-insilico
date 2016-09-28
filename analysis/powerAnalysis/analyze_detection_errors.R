# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

source('powerAnalysis/lib.R')

# Constants

kGroundtruthPath <- 'powerAnalysis/data/simChrom_groundtruth.txt'
kDetectionsPattern <- 'results/mcmc_detections_powerAnalysis_chrom%02d_pm%d.txt'
kReplicates <- 1:10
kPm <- 0:3
kDetectionPowerPattern <- 'powerAnalysis/output_detection_power_pm%d.RData'
kDetectionPowerPlots <- 'powerAnalysis/plots_power_pm%d.pdf'
theme_set(theme_bw())

# Load groundtruth data
groundtruth <- read.csv(kGroundtruthPath)
groundtruth$rep <- groundtruth$rep + 1
groundtruth <- groundtruth[order(groundtruth$rep, groundtruth$pos), ]

# Label primaries with design offsets
groundtruth$offset.primary <- 0
is.na(groundtruth$offset.primary) <- groundtruth$type > 0
ind <- which(groundtruth$type == 0 & groundtruth$eff.magnitude < 1)
groundtruth$offset.primary[ind] <- groundtruth$offset[ind + 1]

for (pm in kPm) {
  detections.list <- list()
  for (r in kReplicates) {
    # Get indices to for relevant subset of groundtruth
    subset.rep <- which(groundtruth$rep == r)
    
    # Load detections for each replicate
    detections.list[[r]] <- read.table(sprintf(kDetectionsPattern, r, pm),
                                       header=TRUE)
    detections.list[[r]]$rep <- r
    
    # Find nearest-neighbor positions
    matched.groundtruth <- MatchPositions(detections.list[[r]]$pos, 
                                          groundtruth$pos[subset.rep])
    matched.groundtruth.errors <- detections.list[[r]]$pos - matched.groundtruth
    
    matched.detections <- MatchPositions(groundtruth$pos[subset.rep],
                                         detections.list[[r]]$pos)
    matched.detections.errors <- groundtruth$pos[subset.rep] - matched.detections
    
    # Build data frames for detailed analysis
    groundtruth$error[subset.rep] <- matched.detections.errors  
    detections.list[[r]]$error <- matched.groundtruth.errors
  }
  
  # Combined detections information
  detections <- Reduce(rbind, detections.list)
  rm(detections.list)
  gc()
  
  # Compute error indicators
  groundtruth$detected.alternative <- (abs(groundtruth$error) <=
                                         groundtruth$offset/2)
  groundtruth$detected.primary <- (abs(groundtruth$error) < 5)
  
  # Compute aggregate error summaries
  
  # Primary
  power.primary <- ddply(groundtruth[groundtruth$type == 0, ],
                         .(eff.magnitude, coverage, offset.primary),
                         summarise,
                         p.detected.primary=mean(detected.primary),
                         mean.err=mean(abs(error)),
                         se.err=sd(abs(error)) / sqrt(length(error)),
                         n=length(detected.primary))
  power.primary.by.eff.magnitude <- ddply(
    power.primary, .(eff.magnitude),summarise,
    p.detected.primary=mean(p.detected.primary),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  power.primary.by.coverage <- ddply(
    power.primary, .(coverage), summarise,
    p.detected.primary=mean(p.detected.primary),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  power.primary.by.offset <- ddply(
    power.primary, .(offset.primary), summarise,
    p.detected.primary=mean(p.detected.primary),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  
  # Alternative
  power.alternative <- ddply(groundtruth[groundtruth$type > 0, ],
                             .(eff.magnitude, coverage, offset),
                             summarise,
                             p.detected.alternative=mean(detected.alternative),
                             mean.err=mean(abs(error)),
                             se.err=sd(abs(error)) / sqrt(length(error)),
                             n=length(detected.alternative))
  power.alternative.by.eff.magnitude <- ddply(
    power.alternative, .(eff.magnitude),summarise,
    p.detected.alternative=mean(p.detected.alternative),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  power.alternative.by.coverage <- ddply(
    power.alternative, .(coverage), summarise,
    p.detected.alternative=mean(p.detected.alternative),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  power.alternative.by.offset <- ddply(
    power.alternative, .(offset), summarise,
    p.detected.alternative=mean(p.detected.alternative),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  
  
  # Regression-based summaries
  
  power.primary$eff.magnitude.adj <- with(power.primary,
                                          eff.magnitude - min(eff.magnitude))
  power.primary$coverage.adj <- with(power.primary,
                                     coverage - min(coverage))
  lm.primary.power <- lm(
    p.detected.primary ~ eff.magnitude.adj * coverage.adj * offset.primary,
    data=power.primary)
  glm.primary.power <- glm(
    p.detected.primary ~ eff.magnitude.adj * coverage.adj * offset.primary,
    data=power.primary,
    family=binomial(link='logit'),
    weights=n)
  
  power.alternative$eff.magnitude.adj <- with(power.alternative,
                                              eff.magnitude - min(eff.magnitude))
  power.alternative$coverage.adj <- with(power.alternative,
                                         coverage - min(coverage))
  
  lm.alternative.power <- lm(
    p.detected.alternative ~ eff.magnitude.adj * coverage.adj * offset,
    data=power.alternative)
  glm.alternative.power <- glm(
    p.detected.alternative ~ eff.magnitude.adj * coverage.adj * offset,
    data=power.alternative,
    family=binomial(link='logit'),
    weights=n)
  
  # Build basic summary plots
  
  # Primaries
  
  plot.list.primary <- list()
  
  plot.list.primary[[1]] <- qplot(x=coverage, y=p.detected.primary,
                                  data=power.primary,
                                  geom=c('line', 'point'), ylim=c(0, 1))
  plot.list.primary[[1]] <- plot.list.primary[[1]] + facet_grid(
    offset.primary ~ eff.magnitude, labeller=function(var, value)
      round(value,3))
  plot.list.primary[[1]] <- plot.list.primary[[1]] + labs(list(
    x='Coverage',
    y='Detection probability',
    title='Primary position detection probability by offset and effective magnitude'))
  
  plot.list.primary[[2]] <- qplot(x=coverage, y=p.detected.primary,
                                  data=power.primary.by.coverage,
                                  geom=c('line', 'point'), ylim=c(0, 1))
  plot.list.primary[[2]] <- plot.list.primary[[2]] + labs(list(
    x='Coverage',
    y='Detection probability',
    title='Primary position detection probability vs. coverage'))
  
  plot.list.primary[[3]] <- qplot(x=eff.magnitude, y=p.detected.primary,
                                  data=power.primary.by.eff.magnitude,
                                  geom=c('line', 'point'), ylim=c(0, 1))
  plot.list.primary[[3]] <- plot.list.primary[[3]] + labs(list(
    x='Effective magnitude',
    y='Detection probability',
    title='Primary position detection probability vs. effective magnitude'))
  
  plot.list.primary[[4]] <- qplot(x=offset.primary, y=p.detected.primary,
                                  data=power.primary.by.offset,
                                  geom=c('line', 'point'), ylim=c(0, 1))
  plot.list.primary[[4]] <- plot.list.primary[[4]] + labs(list(
    x='Offset',
    y='Detection probability',
    title='Primary position detection probability vs. offset'))
  
  
  plot.list.primary[[5]] <- qplot(
    x=coverage, y=mean.err,
    ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
    data=power.primary,
    geom=c('line', 'errorbar', 'point'))
  plot.list.primary[[5]] <- plot.list.primary[[5]] + facet_grid(
    offset.primary ~ eff.magnitude, labeller=function(var, value)
      round(value,3))
  plot.list.primary[[5]]<- plot.list.primary[[5]] + labs(list(
    x='Coverage',
    y='Mean absolute error',
    title='Primary position errors by offset and effective magnitude'))
  
  plot.list.primary[[6]] <- qplot(
    x=coverage, y=mean.err,
    ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
    data=power.primary.by.coverage,
    geom=c('line', 'errorbar', 'point'))
  plot.list.primary[[6]] <- plot.list.primary[[6]] + labs(list(
    x='Coverage',
    y='Mean absolute error',
    title='Primary position position errors vs. coverage'))
  
  plot.list.primary[[7]] <- qplot(
    x=eff.magnitude, y=mean.err,
    ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
    data=power.primary.by.eff.magnitude,
    geom=c('line', 'errorbar', 'point'))
  plot.list.primary[[7]] <- plot.list.primary[[7]] + labs(list(
    x='Effective magnitude',
    y='Mean absolute error',
    title='Primary position position errors vs. effective magnitude'))
  
  plot.list.primary[[8]] <- qplot(
    x=offset.primary, y=mean.err,
    ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
    data=power.primary.by.offset,
    geom=c('line', 'errorbar', 'point'))
  plot.list.primary[[8]] <- plot.list.primary[[8]] + labs(list(
    x='Offset',
    y='Mean absolute error',
    title='Primary position position errors vs. effective magnitude'))
  
  # Alternative
  
  plot.list.alternative <- list()
  
  plot.list.alternative[[1]] <- qplot(x=coverage, y=p.detected.alternative,
                                      data=power.alternative,
                                      geom=c('line', 'point'), ylim=c(0, 1))
  plot.list.alternative[[1]] <- plot.list.alternative[[1]] + facet_grid(
    offset ~ eff.magnitude, labeller=function(var, value) round(value,3))
  plot.list.alternative[[1]] <- plot.list.alternative[[1]] + labs(list(
    x='Coverage',
    y='Detection probability',
    title='Alternative position detection probability by offset and effective magnitude'))
  
  plot.list.alternative[[2]] <- qplot(x=coverage, y=p.detected.alternative,
                                      data=power.alternative.by.coverage,
                                      geom=c('line', 'point'), ylim=c(0, 1))
  plot.list.alternative[[2]] <- plot.list.alternative[[2]] + labs(list(
    x='Coverage',
    y='Detection probability',
    title='Alternative position detection probability vs. coverage'))
  
  plot.list.alternative[[3]] <- qplot(x=eff.magnitude, y=p.detected.alternative,
                                      data=power.alternative.by.eff.magnitude,
                                      geom=c('line', 'point'), ylim=c(0, 1))
  plot.list.alternative[[3]] <- plot.list.alternative[[3]] + labs(list(
    x='Effective magnitude',
    y='Detection probability',
    title='Alternative position detection probability vs. effective magnitude'))
  
  plot.list.alternative[[4]] <- qplot(x=offset, y=p.detected.alternative,
                                      data=power.alternative.by.offset,
                                      geom=c('line', 'point'), ylim=c(0, 1))
  plot.list.alternative[[4]] <- plot.list.alternative[[4]] + labs(list(
    x='Offset',
    y='Detection probability',
    title='Alternative position detection probability vs. offset'))
  
  
  plot.list.alternative[[5]] <- qplot(
    x=coverage, y=mean.err,
    ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
    data=power.alternative,
    geom=c('line', 'errorbar', 'point'))
  plot.list.alternative[[5]] <- plot.list.alternative[[5]] + facet_grid(
    offset ~ eff.magnitude, labeller=function(var, value) round(value,3))
  plot.list.alternative[[5]]<- plot.list.alternative[[5]] + labs(list(
    x='Coverage',
    y='Mean absolute error',
    title='Alternative position errors by offset and effective magnitude'))
  
  plot.list.alternative[[6]] <- qplot(
    x=coverage, y=mean.err,
    ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
    data=power.alternative.by.coverage,
    geom=c('line', 'errorbar', 'point'))
  plot.list.alternative[[6]] <- plot.list.alternative[[6]] + labs(list(
    x='Coverage',
    y='Mean absolute error',
    title='Alternative position position errors vs. coverage'))
  
  plot.list.alternative[[7]] <- qplot(
    x=eff.magnitude, y=mean.err,
    ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
    data=power.alternative.by.eff.magnitude,
    geom=c('line', 'errorbar', 'point'))
  plot.list.alternative[[7]] <- plot.list.alternative[[7]] + labs(list(
    x='Effective magnitude',
    y='Mean absolute error',
    title='Alternative position position errors vs. effective magnitude'))
  
  plot.list.alternative[[8]] <- qplot(
    x=offset, y=mean.err,
    ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
    data=power.alternative.by.offset,
    geom=c('line', 'errorbar', 'point'))
  plot.list.alternative[[8]] <- plot.list.alternative[[8]] + labs(list(
    x='Offset',
    y='Mean absolute error',
    title='Alternative position position errors vs. offset'))
  
  # Save plots
  pdf(sprintf(kDetectionPowerPlots, pm), 11, 11*0.62, onefile=TRUE)
  l_ply(plot.list.primary, print)
  l_ply(plot.list.alternative, print)
  dev.off()
  
  # Save results
  save(groundtruth, detections, power.primary, power.alternative,
       file=sprintf(kDetectionPowerPattern, pm))
}
