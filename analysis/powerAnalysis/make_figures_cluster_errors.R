# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(grid)
library(gridExtra)

source("powerAnalysis/lib.R")

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


# Load cluster results
cluster.env <- new.env()
load("powerAnalysis/output_cluster_power.RData", envir=cluster.env)

# Load Parzen window results
parzen.env <- new.env()
load("powerAnalysis/output_parzen_window_power.RData", envir=parzen.env)

# Load NOrMAL results
normal.env <- new.env()
normal.env$power.cluster <- read.table("comparisons/out/NOrMAL.tsv",
                                       header=TRUE, sep="\t")
normal.env$power.cluster$coverage <-
  as.numeric(as.character(normal.env$power.cluster$coverage))

# Combine results into single data.frame
cluster.env$power.cluster$Method <- "Cluster estimand"
parzen.env$power.cluster$Method <- "Parzen window"
normal.env$power.cluster$Method <- "NOrMAL"
p.detected.df <- rbind(cluster.env$power.cluster, parzen.env$power.cluster,
                       normal.env$power.cluster)
p.detected.df$coverage <- round(p.detected.df$coverage, 2)
p.detected.df$se.power <- with(
  p.detected.df, sqrt(p.detected.primary * (1 - p.detected.primary) / n))
power.by.eff.magnitude <- ddply(
  p.detected.df, .(Method, eff.magnitude), 
  summarise,
  power=mean(p.detected.primary),
  se.power=sqrt(mean(se.power^2) / length(p.detected.primary)),
  mean.err=mean(mean.err),
  se.err=sqrt(mean(se.err^2) / length(se.err)))
power.by.offset <- ddply(
  p.detected.df, .(Method, offset.primary), 
  summarise,
  power=mean(p.detected.primary),
  se.power=sqrt(mean(se.power^2) / length(p.detected.primary)),
  mean.err=mean(mean.err),
  se.err=sqrt(mean(se.err^2) / length(se.err)))
power.by.coverage <- ddply(
  p.detected.df, .(Method, coverage), 
  summarise,
  power=mean(p.detected.primary),
  se.power=sqrt(mean(se.power^2) / length(p.detected.primary)),
  mean.err=mean(mean.err),
  se.err=sqrt(mean(se.err^2) / length(se.err)))

# Build plots
blank_legend <- theme(legend.key=element_rect(color="white"),
                      legend.text=element_blank(),
                      legend.title=element_blank())

panel.list <- list()

panel.list[[1]] <- qplot(
  x=eff.magnitude, y=power,
  ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
  data=power.by.eff.magnitude,
  color=Method, shape=Method,
  geom=c("line", "point", "errorbar"), ylim=c(0, 1))
panel.list[[1]] <- panel.list[[1]] + labs(list(
  x="Effective magnitude",
  y="Power"))


panel.list[[2]] <- qplot(
  x=offset.primary, y=power,
  ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
  data=power.by.offset,
  color=Method, shape=Method,
  geom=c("line", "point", "errorbar"), ylim=c(0, 1))
panel.list[[2]] <- panel.list[[2]] + labs(list(
  x="Offset",
  y="Power"))

panel.list[[3]] <- qplot(
  x=coverage, y=power,
  ymin=power - 1.96 * se.power, ymax=power + 1.96 * se.power,
  data=power.by.coverage,
  color=Method, shape=Method,
  geom=c("line", "point", "errorbar"), ylim=c(0, 1))
panel.list[[3]] <- panel.list[[3]] + labs(list(
  x="Coverage",
  y="Power"))

legend.v <- GetLegend(panel.list[[1]])
legend.h <- GetLegend(panel.list[[1]] + theme(legend.direction="horizontal"))

theme_set(theme_bw())

pdf("powerAnalysis/figure_power_cluster_2-panel.pdf", 10, 3.1)

grid.arrange(panel.list[[1]] + theme(legend.position="none",
                                     axis.title.y=element_blank()),
             panel.list[[2]] + theme(legend.position="none",
                                     axis.title.y=element_blank()),
             legend.v,
             left="Power",
             widths=unit.c(unit(0.5, "npc") - 0.5*legend.v$width[2],
                           unit(0.5, "npc") - 0.5*legend.v$width[2],
                           legend.v$width[2]),
             ncol=3)

dev.off()

pdf("powerAnalysis/figure_power_cluster_3-panel.pdf", 10, 4)

grid.arrange(arrangeGrob(panel.list[[1]] + theme(legend.position="none",
                                                 axis.title.y=element_blank()),
                         panel.list[[2]] + theme(legend.position="none",
                                                 axis.title.y=element_blank()),
                         panel.list[[3]] + theme(legend.position="none",
                                                 axis.title.y=element_blank()),
                         ncol=3, left="Power"),
             legend.h,
             widths=unit.c(unit(1, "npc")),
             heights=unit.c(unit(1, "npc") - legend.h$heights[2],
                            legend.h$heights[2]),
             nrow=2, ncol=1)

dev.off()
