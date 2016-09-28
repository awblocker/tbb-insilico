# Load libraries
library(plyr)
library(reshape2)
library(stringr)
library(grid)
library(gridExtra)
library(Rcpp)
library(ggplot2)
library(ascii)

source('powerAnalysis/lib.R')

# Set constants
kInputPaths <- list(
  groundtruth='powerAnalysis/output_cluster_summaries_groundtruth.RData',
  Model='powerAnalysis/output_cluster_summaries_model.RData',
  Read='powerAnalysis/output_cluster_summaries_parzen.RData')
kPlotPathPattern <- 'powerAnalysis/figure_cluster_summaries.%s'
kVars <- c('localization', 'sparsity', 'structure')
kTextPath <- 'powerAnalysis/results_cluster_summaries.txt'
kOutputPath <- 'powerAnalysis/output_cluster_summary_analysis.RData'


# Function definitions


MatchToGroundtruth <- function(clusters, groundtruth) {
  clusters.list <- list()
  
  clusters$matched.center <- 0
  clusters$matched.localization <- 0
  clusters$matched.sparsity <- 0
  clusters$matched.structure <- 0
  
  for (r in unique(clusters$rep)) {
    # Get indices to for relevant subset of groundtruth
    subset.groundtruth.rep <- which(groundtruth$rep == r)
    subset.clusters.rep <- which(clusters$rep == r)
    
    # Find nearest-neighbor positions
    matched.groundtruth <- MatchPositions(
      clusters$center[subset.clusters.rep], 
      groundtruth$center[subset.groundtruth.rep])
    matched.indices <- match(matched.groundtruth,
                             groundtruth$center[subset.groundtruth.rep])
    
    # Add matched information to clusters
    clusters$matched.center[subset.clusters.rep] <- matched.groundtruth
    if ('localization' %in% names(clusters))
      clusters$matched.localization[subset.clusters.rep] <- groundtruth[[
        'localization']][subset.groundtruth.rep[matched.indices]]
    if ('sparsity' %in% names(clusters))
      clusters$matched.sparsity[subset.clusters.rep] <- groundtruth[[
        'sparsity']][subset.groundtruth.rep[matched.indices]]
    if ('structure' %in% names(clusters))
      clusters$matched.structure[subset.clusters.rep] <- groundtruth[[
        'structure']][subset.groundtruth.rep[matched.indices]]
  }
  
  return(clusters)
}


# Script


# Load data into environments
data.env.list <- lapply(kInputPaths, function(path) {
  io.env <- new.env(hash=TRUE)
  load(path, envir=io.env)
  return(io.env)
})

# Run matching
results.list <- lapply(
  names(data.env.list)[names(data.env.list) != 'groundtruth'],
  function(data.name) {
    out <- MatchToGroundtruth(data.env.list[[data.name]]$clusters,
                              data.env.list$groundtruth$clusters)
    out$Method <- data.name
    return(out)
  })
names(results.list) <- names(data.env.list)[
  names(data.env.list) != 'groundtruth']

common.cols <- Reduce(intersect, lapply(results.list, names))
for (ii in 1:length(results.list))
  results.list[[ii]] <- results.list[[ii]][, common.cols]

results <- Reduce(rbind, results.list)

# Compute aggregate summaries

results.long <- reshape(results, varying=list(kVars, str_c('matched.', kVars)),
                        v.names=c('value', 'matched'),
                        idvar=c('rep', 'center', 'Method'),
                        timevar='var', times=kVars,
                        direction='long')

summaries.regression <- ddply(
  results.long, .(var, Method), function(dd) {
    model <- lm(value ~ matched, data=dd)
    c(coef(model), r.squared=summary(model)$r.squared)
  })

# Build plots

plots.by.method <- dlply(results, .(Method), function(dd) {
  plt.structure <- qplot(
    x=matched.structure, y=structure, data=dd,
    geom='point', size=I(1), alpha=I(0.2),
    ylim=c(0, 1),
    xlab='Structure', ylab='Estimated structure')
  plt.structure <- plt.structure + geom_abline(
    intercept=0, slope=1, colour=I('red'))
  
  plt.localization <- qplot(
    x=matched.localization, y=localization, data=dd,
    geom='point', size=I(1), alpha=I(0.2),
    ylim=c(-1, 1),
    xlab='Localization', ylab='Estimated localization')
  plt.localization <- plt.localization + geom_abline(
    intercept=0, slope=1, colour=I('red'))
  
  plt.sparsity <- qplot(
    x=matched.sparsity, y=sparsity, data=dd,
    ylim=c(0, 1),
    geom='point', size=I(1), alpha=I(0.2),
    xlab='Sparsity', ylab='Estimated sparsity',
    position='jitter', w=0.01, h=0)
  plt.sparsity <- plt.sparsity + geom_abline(
    intercept=0, slope=1, colour=I('red'))
  
  return(list(Localization=plt.localization,
              Sparsity=plt.sparsity,
              Structure=plt.structure))
})

# Assemble and output combined plots

theme_set(theme_bw())

tiff(sprintf(kPlotPathPattern, 'tiff'),
     compression='lzw',
     width=10, height=4, units='in', res=300)
grid.arrange(
  arrangeGrob(
    plots.by.method$Model$Localization +
      ggtitle('Localization') + theme(axis.title.x=element_blank()),
    plots.by.method$Model$Sparsity +
      ggtitle('Sparsity') + theme(axis.title.x=element_blank()),
    plots.by.method$Model$Structure +
      ggtitle('Structure') + theme(axis.title.x=element_blank()),
    ncol=3, left='Model-based'),
  arrangeGrob(
    plots.by.method$Read$Localization,
    plots.by.method$Read$Sparsity,
    plots.by.method$Read$Structure,
    ncol=3, left='Read-based'),
  widths=unit.c(unit(1, "npc")),
  nrow=2, ncol=1)
dev.off()
system(sprintf('convert -density 300 %s %s',
               sprintf(kPlotPathPattern, 'tiff'),
               sprintf(kPlotPathPattern, 'pdf')))

# Save text output
print(ascii(summaries.regression,
            include.rownames=FALSE, include.colnames=TRUE,
            header=FALSE, format='g', digits=4),
      format='rest',
      file=kTextPath)
