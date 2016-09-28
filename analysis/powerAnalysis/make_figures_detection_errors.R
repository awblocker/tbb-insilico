# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(grid)
library(gridExtra)

source('powerAnalysis/lib.R')

# Set constants
kPm <- 0:3
kDetectionPowerPattern <- 'powerAnalysis/output_detection_power_pm%d.RData'
kPlotPathPattern <- 'powerAnalysis/figure_power_%s-pm%d%s.pdf'

# Function definitions
VpLayout <- function(x, y) {
  viewport(layout.pos.row=x, layout.pos.col=y)
}

GetLegend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# Script

for (pm in kPm) {
  # Load detection results
  load(sprintf(kDetectionPowerPattern, pm))
  

  # Process primary power statistics ----------------------------------------
  power.primary$se.power <- with(
    power.primary, sqrt(p.detected.primary * (1 - p.detected.primary) / n))
  power.primary.by.eff.magnitude <- ddply(
    power.primary, .(eff.magnitude), 
    summarise,
    power=mean(p.detected.primary),
    se.power=sqrt(mean(se.power^2) / length(p.detected.primary)),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  power.primary.by.offset <- ddply(
    power.primary, .(offset.primary), 
    summarise,
    power=mean(p.detected.primary),
    se.power=sqrt(mean(se.power^2) / length(p.detected.primary)),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  power.primary.by.coverage <- ddply(
    power.primary, .(coverage), 
    summarise,
    power=mean(p.detected.primary),
    se.power=sqrt(mean(se.power^2) / length(p.detected.primary)),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  
  

  # Process alternative power statistics ------------------------------------
  power.alternative$se.power <- with(
    power.alternative, sqrt(p.detected.alternative *
                              (1 - p.detected.alternative) / n))
  power.alternative.by.eff.magnitude <- ddply(
    power.alternative, .(eff.magnitude), 
    summarise,
    power=mean(p.detected.alternative),
    se.power=sqrt(mean(se.power^2) / length(p.detected.alternative)),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  power.alternative.by.offset <- ddply(
    power.alternative, .(offset), 
    summarise,
    power=mean(p.detected.alternative),
    se.power=sqrt(mean(se.power^2) / length(p.detected.alternative)),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))
  power.alternative.by.coverage <- ddply(
    power.alternative, .(coverage), 
    summarise,
    power=mean(p.detected.alternative),
    se.power=sqrt(mean(se.power^2) / length(p.detected.alternative)),
    mean.err=mean(mean.err),
    se.err=sqrt(mean(se.err^2) / length(se.err)))

  # Build primary plots
  panel.primary.list <- list()
  
  panel.primary.list[[1]] <- qplot(
    x=eff.magnitude, y=power,
    ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
    data=power.primary.by.eff.magnitude,
    geom=c('line', 'point', 'errorbar'), ylim=c(0, 1))
  panel.primary.list[[1]] <- panel.primary.list[[1]] + labs(list(
    x='Effective magnitude',
    y='Power'))
  
  
  panel.primary.list[[2]] <- qplot(
    x=offset.primary, y=power,
    ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
    data=power.primary.by.offset,
    geom=c('line', 'point', 'errorbar'), ylim=c(0, 1))
  panel.primary.list[[2]] <- panel.primary.list[[2]] + labs(list(
    x='Offset',
    y='Power'))
  
  panel.primary.list[[3]] <- qplot(
    x=coverage, y=power,
    ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
    data=power.primary.by.coverage,
    geom=c('line', 'point', 'errorbar'), ylim=c(0, 1))
  panel.primary.list[[3]] <- panel.primary.list[[3]] + labs(list(
    x='Coverage',
    y='Power'))  
  
  # Output primary plots
  
  theme_set(theme_bw())
  
  pdf(sprintf(kPlotPathPattern, 'primary', pm, '_2-panel'), 10, 3.1)
  
  grid.arrange(panel.primary.list[[1]] + theme(legend.position="none",
                                       axis.title.y=element_blank()),
               panel.primary.list[[2]] + theme(legend.position="none",
                                       axis.title.y=element_blank()),
               left='Power',
               widths=unit.c(unit(0.5, "npc"),  unit(0.5, "npc")),
               ncol=2
               )
  
  dev.off()
  
  pdf(sprintf(kPlotPathPattern, 'primary', pm, '_3-panel'), 10, 3.1)
  
  grid.arrange(arrangeGrob(
    panel.primary.list[[1]] + theme(legend.position="none",
                                    axis.title.y=element_blank()),
    panel.primary.list[[2]] + theme(legend.position="none",
                                    axis.title.y=element_blank()),
    panel.primary.list[[3]] + theme(legend.position="none",
                                    axis.title.y=element_blank()),
    ncol=3, left='Power'),
               widths=unit.c(unit(1, "npc")),
               nrow=1, ncol=1)
  
  dev.off()
  
  # Build alternative plots
  
  panel.alternative.list <- list()
  
  panel.alternative.list[[1]] <- qplot(
    x=eff.magnitude, y=power,
    ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
    data=power.alternative.by.eff.magnitude,
    geom=c('line', 'point', 'errorbar'), ylim=c(0, 1))
  panel.alternative.list[[1]] <- panel.alternative.list[[1]] + labs(list(
    x='Effective magnitude',
    y='Power'))
  
  
  panel.alternative.list[[2]] <- qplot(
    x=offset, y=power,
    ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
    data=power.alternative.by.offset,
    geom=c('line', 'point', 'errorbar'), ylim=c(0, 1))
  panel.alternative.list[[2]] <- panel.alternative.list[[2]] + labs(list(
    x='Offset',
    y='Power'))
  
  panel.alternative.list[[3]] <- qplot(
    x=coverage, y=power,
    ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
    data=power.alternative.by.coverage,
    geom=c('line', 'point', 'errorbar'), ylim=c(0, 1))
  panel.alternative.list[[3]] <- panel.alternative.list[[3]] + labs(list(
    x='Coverage',
    y='Power'))  
  
  # Output alternative plots
  
  theme_set(theme_bw())
  
  pdf(sprintf(kPlotPathPattern, 'alternative', pm, '_2-panel'), 10, 3.1)
  
  grid.arrange(panel.alternative.list[[1]] + theme(legend.position="none",
                                       axis.title.y=element_blank()),
               panel.alternative.list[[2]] + theme(legend.position="none",
                                       axis.title.y=element_blank()),
               left='Power',
               widths=unit.c(unit(0.5, "npc"), unit(0.5, "npc")),
               ncol=2
  )
  
  dev.off()
  
  pdf(sprintf(kPlotPathPattern, 'alternative', pm, '_3-panel'), 10, 3.1)
  
  grid.arrange(arrangeGrob(
    panel.alternative.list[[1]] + theme(legend.position="none",
                                        axis.title.y=element_blank()),
    panel.alternative.list[[2]] + theme(legend.position="none",
                                        axis.title.y=element_blank()),
    panel.alternative.list[[3]] + theme(legend.position="none",
                                        axis.title.y=element_blank()),
    ncol=3, left='Power'),
               widths=unit.c(unit(1, "npc")),
               nrow=1, ncol=1)
  
  dev.off()
  
  # Combined plot with primary and alternative
  
  pdf(sprintf(kPlotPathPattern, 'combined', pm, '_3-panel'), 10, 4)
  
  grid.arrange(
    arrangeGrob(
      panel.primary.list[[1]] + theme(legend.position="none",
                                      axis.title.y=element_blank()),
      panel.primary.list[[2]] + theme(legend.position="none",
                                      axis.title.y=element_blank()),
      panel.primary.list[[3]] + theme(legend.position="none",
                                      axis.title.y=element_blank()),
      ncol=3, left='Primary'),
    arrangeGrob(
      panel.alternative.list[[1]] + theme(legend.position="none",
                                      axis.title.y=element_blank()),
      panel.alternative.list[[2]] + theme(legend.position="none",
                                      axis.title.y=element_blank()),
      panel.alternative.list[[3]] + theme(legend.position="none",
                                      axis.title.y=element_blank()),
      ncol=3, left='Alternative'),
    left='Power', widths=unit.c(unit(1, "npc")),
    nrow=2, ncol=1)
  
  dev.off()
}
