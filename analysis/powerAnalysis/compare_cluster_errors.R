# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

source('powerAnalysis/lib.R')

# Load cluster results
cluster.env <- new.env()
load('powerAnalysis/output_cluster_power.RData', envir=cluster.env)

# Load Parzen window results
parzen.env <- new.env()
load('powerAnalysis/output_parzen_window_power.RData', envir=parzen.env)

# Combine results into single data.frame
cluster.env$power.cluster$method <- 'Model-based'
parzen.env$power.cluster$method <- 'Parzen window'
power.cluster <- rbind(cluster.env$power.cluster, parzen.env$power.cluster)
power.by.eff.magnitude <- ddply(power.cluster, .(method, eff.magnitude), 
                                summarise,
                                p.detected.primary=mean(p.detected.primary),
                                mean.err=mean(mean.err),
                                se.err=sqrt(mean(se.err^2) / length(se.err)))
power.by.coverage <- ddply(power.cluster, .(method, coverage), 
                           summarise,
                           p.detected.primary=mean(p.detected.primary),
                           mean.err=mean(mean.err),
                           se.err=sqrt(mean(se.err^2) / length(se.err)))
power.by.offset <- ddply(power.cluster, .(method, offset.primary), 
                           summarise,
                           p.detected.primary=mean(p.detected.primary),
                           mean.err=mean(mean.err),
                           se.err=sqrt(mean(se.err^2) / length(se.err)))

# Build plots
plot.list <- list()

plot.list[[1]] <- qplot(x=coverage, y=p.detected.primary,
                        data=power.cluster, color=method,
                        geom=c('line', 'point'), ylim=c(0, 1))
plot.list[[1]] <- plot.list[[1]] + facet_grid(
  offset.primary ~ eff.magnitude, labeller=function(var, value) round(value,3))
plot.list[[1]] <- plot.list[[1]] + labs(list(
  x='Coverage',
  y='Detection probability',
  title='Model-based and Parzen window primary position detection probability\nby offset and effective magnitude'))

plot.list[[2]] <- qplot(x=coverage, y=p.detected.primary,
                        data=power.by.coverage, color=method,
                        geom=c('line', 'point'), ylim=c(0, 1))
plot.list[[2]] <- plot.list[[2]] + labs(list(
  x='Coverage',
  y='Detection probability',
  title='Model-based and Parzen window primary position detection probability vs. coverage'))

plot.list[[3]] <- qplot(x=eff.magnitude, y=p.detected.primary,
                        data=power.by.eff.magnitude, color=method,
                        geom=c('line', 'point'), ylim=c(0, 1))
plot.list[[3]] <- plot.list[[3]] + labs(list(
  x='Effective magnitude',
  y='Detection probability',
  title='Model-based and Parzen window primary position detection probability vs. effective magnitude'))

plot.list[[4]] <- qplot(x=offset.primary, y=p.detected.primary,
                        data=power.by.offset, color=method,
                        geom=c('line', 'point'), ylim=c(0, 1))
plot.list[[4]] <- plot.list[[4]] + labs(list(
  x='Offset',
  y='Detection probability',
  title='Model-based and Parzen window primary position detection probability vs. offset'))


plot.list[[5]] <- qplot(x=coverage, y=mean.err,
                        ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
                        data=power.cluster, color=method,
                        geom=c('line', 'errorbar', 'point'))
plot.list[[5]] <- plot.list[[5]] + facet_grid(
  offset.primary ~ eff.magnitude, labeller=function(var, value) round(value,3))
plot.list[[5]]<- plot.list[[5]] + labs(list(
  x='Coverage',
  y='Mean absolute error',
  title='Model-based and Parzen window primary position errors\nby offset and effective magnitude'))

plot.list[[6]] <- qplot(x=coverage, y=mean.err,
                        ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
                        data=power.by.coverage, color=method,
                        geom=c('line', 'errorbar', 'point'))
plot.list[[6]] <- plot.list[[6]] + labs(list(
  x='Coverage',
  y='Mean absolute error',
  title='Model-based and Parzen window primary position position errors vs. coverage'))

plot.list[[7]] <- qplot(x=eff.magnitude, y=mean.err,
                        ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
                        data=power.by.eff.magnitude, color=method,
                        geom=c('line', 'errorbar', 'point'))
plot.list[[7]] <- plot.list[[7]] + labs(list(
  x='Effective magnitude',
  y='Mean absolute error',
  title='Model-based and Parzen window primary position position errors vs. effective magnitude'))

plot.list[[8]] <- qplot(x=offset.primary, y=mean.err,
                        ymin=mean.err - 2 * se.err, ymax=mean.err + 2 * se.err,
                        data=power.by.offset, color=method,
                        geom=c('line', 'errorbar', 'point'))
plot.list[[8]] <- plot.list[[8]] + labs(list(
  x='Offset',
  y='Mean absolute error',
  title='Model-based and Parzen window primary position position errors vs. effective magnitude'))

theme_set(theme_bw())
pdf('powerAnalysis/plots_compare_power.pdf', 11, 11*0.62, onefile=TRUE)
l_ply(plot.list, print)
dev.off()
