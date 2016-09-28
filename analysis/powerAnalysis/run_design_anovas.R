# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

source('powerAnalysis/lib.R')

# Set constants
kInputs <- list(cluster='powerAnalysis/output_cluster_power.RData',
                normal='comparisons/out/NOrMAL.Rdata',
                parzen='powerAnalysis/output_parzen_window_power.RData',
                local='powerAnalysis/output_detection_power_pm3.RData')
kOutputPath <- 'powerAnalysis/results_design_anovas.RData'
kGeneLength <- 3501

# Function definitions

LoadIntoNewEnv <- function(path) {
  env.data <- new.env(hash=TRUE)
  load(path, envir=env.data)
  return(env.data)
}


# Load results
results.list <- lapply(kInputs, LoadIntoNewEnv)

# Extract data for ANOVAs
data.list <- list(
  cluster=results.list$cluster$groundtruth[
    results.list$cluster$groundtruth$type == 0, ],
  normal=results.list$normal$groundtruth[
    results.list$normal$groundtruth$type == 0, ],
  parzen=results.list$parzen$groundtruth[
    results.list$parzen$groundtruth$type == 0, ])
data.list$primary <- results.list$local$groundtruth[
  results.list$local$groundtruth$type == 0, ]
data.list$alternative <- results.list$local$groundtruth[
  results.list$local$groundtruth$type > 0, ]

# Add gene identifiers to data
for (ii in 1:length(data.list)) {
  data.list[[ii]]$gene <- as.integer((data.list[[ii]]$pos - 1) %/% kGeneLength)
}

# Aggregate mean errors by gene
data.list$cluster <- ddply(
  data.list$cluster, .(rep, gene), summarise,
  offset.primary=offset.primary[1],
  coverage=coverage[1],
  eff.magnitude=eff.magnitude[1],
  magnitude=magnitude[1],
  p.detected.primary=mean(detected.primary),
  mean.abs.error=mean(abs(error)),
  se.abs.error=sd(abs(error)) / sqrt(length(error)),
  n=length(error))
data.list$normal <- ddply(
  data.list$normal, .(rep, gene), summarise,
  offset.primary=as.numeric(as.character(offset.primary[1])),
  coverage=as.numeric(as.character(coverage[1])),
  eff.magnitude=as.numeric(as.character(eff.magnitude[1])),
  magnitude=as.numeric(as.character(magnitude[1])),
  p.detected.primary=mean(detected.primary),
  mean.abs.error=mean(abs(error)),
  se.abs.error=sd(abs(error)) / sqrt(length(error)),
  n=length(error))
data.list$parzen <- ddply(
  data.list$parzen, .(rep, gene), summarise,
  offset.primary=offset.primary[1],
  coverage=coverage[1],
  eff.magnitude=eff.magnitude[1],
  magnitude=magnitude[1],
  p.detected.primary=mean(detected.primary),
  mean.abs.error=mean(abs(error)),
  se.abs.error=sd(abs(error)) / sqrt(length(error)),
  n=length(error))
data.list$primary <- ddply(
  data.list$primary, .(rep, gene), summarise,
  offset.primary=offset.primary[1],
  coverage=coverage[1],
  eff.magnitude=eff.magnitude[1],
  magnitude=magnitude[1],
  p.detected.primary=mean(detected.primary),
  mean.abs.error=mean(abs(error)),
  se.abs.error=sd(abs(error)) / sqrt(length(error)),
  n=length(error))
data.list$alternative <- ddply(
  data.list$alternative, .(rep, gene), summarise,
  offset=offset[1],
  coverage=coverage[1],
  eff.magnitude=eff.magnitude[1],
  magnitude=magnitude[1],
  p.detected.alternative=mean(detected.alternative),
  mean.abs.error=mean(abs(error)),
  se.abs.error=sd(abs(error)) / sqrt(length(error)),
  n=length(error))

# Estimate inverse-variance weights
w.cluster <- ddply(data.list$cluster,
                  .(offset.primary, coverage, eff.magnitude),
                  summarise, w.hat=1 / mean(se.abs.error^2))
data.list$cluster <- join(data.list$cluster, w.cluster, type='left')
w.normal <- ddply(data.list$normal,
                  .(offset.primary, coverage, eff.magnitude),
                  summarise, w.hat=1 / mean(se.abs.error^2))
data.list$normal <- join(data.list$normal, w.normal, type='left')
w.parzen <- ddply(data.list$parzen,
                  .(offset.primary, coverage, eff.magnitude),
                  summarise, w.hat=1 / mean(se.abs.error^2))
data.list$parzen <- join(data.list$parzen, w.parzen, type='left')
w.primary <- ddply(data.list$primary,
                  .(offset.primary, coverage, eff.magnitude),
                  summarise, w.hat=1 / mean(se.abs.error^2))
data.list$primary <- join(data.list$primary, w.primary, type='left')
w.alternative <- ddply(data.list$alternative,
                   .(offset, coverage, eff.magnitude),
                   summarise, w.hat=1 / mean(se.abs.error^2))
data.list$alternative <- join(data.list$alternative, w.alternative, type='left')

# Build linear models
lm.cluster <- lm(
  mean.abs.error ~ factor(coverage)*factor(offset.primary)*factor(eff.magnitude),
  weights=w.hat,
  data=data.list$cluster, model=FALSE)
lm.normal <- lm(
  mean.abs.error ~ factor(coverage)*factor(offset.primary)*factor(eff.magnitude),
  weights=w.hat,
  data=data.list$normal, model=FALSE)
lm.parzen <- lm(
  mean.abs.error ~ factor(coverage)*factor(offset.primary)*factor(eff.magnitude),
  weights=w.hat,
  data=data.list$parzen, model=FALSE)
lm.primary <- lm(
  mean.abs.error ~ factor(coverage)*factor(offset.primary)*factor(eff.magnitude),
  weights=w.hat,
  data=data.list$primary, model=FALSE)
lm.alternative <- lm(
  mean.abs.error ~ factor(coverage)*factor(offset)*factor(eff.magnitude),
  weights=w.hat,
  data=data.list$alternative, model=FALSE)

# Build GLMs
glm.cluster <- glm(
  p.detected.primary ~ factor(coverage) * factor(offset.primary) *
    factor(eff.magnitude),
  family=binomial(link='logit'),
  weights=n,
  data=results.list$cluster$power.cluster)
glm.normal <- glm(
  p.detected.primary ~ factor(coverage) * factor(offset.primary) *
    factor(eff.magnitude),
  family=binomial(link='logit'),
  weights=n,
  data=results.list$normal$power.cluster)
glm.parzen <- glm(
  p.detected.primary ~ factor(coverage) * factor(offset.primary) *
    factor(eff.magnitude),
  family=binomial(link='logit'),
  weights=n,
  data=results.list$parzen$power.cluster)
glm.primary <- glm(
  p.detected.primary ~ factor(coverage) * factor(offset.primary) *
    factor(eff.magnitude),
  family=binomial(link='logit'),
  weights=n,
  data=results.list$local$power.primary)
glm.alternative <- glm(
  p.detected.alternative ~ factor(coverage) * factor(offset) *
    factor(eff.magnitude),
  family=binomial(link='logit'),
  weights=n,
  data=results.list$local$power.alternative)

# Run ANOVAs
anova.cluster.error <- anova(lm.cluster)
anova.cluster.power <- anova(glm.cluster)

anova.normal.error <- anova(lm.normal)
anova.normal.power <- anova(glm.normal)

anova.parzen.error <- anova(lm.parzen)
anova.parzen.power <- anova(glm.parzen)

anova.primary.error <- anova(lm.primary)
anova.primary.power <- anova(glm.primary)

anova.alternative.error <- anova(lm.alternative)
anova.alternative.power <- anova(glm.alternative)

# Save ANOVAs
save(glm.cluster, anova.cluster.error, anova.cluster.power,
     glm.normal, anova.normal.error, anova.normal.power,
     glm.parzen, anova.parzen.error, anova.parzen.power,
     glm.primary, anova.primary.error, anova.primary.power,
     glm.alternative, anova.alternative.error, anova.alternative.power,
     file=kOutputPath)
