# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

source('powerAnalysis/lib.R')

# Set constants
kInputs <- list(cluster='powerAnalysis/output_cluster_power.RData',
                parzen='powerAnalysis/output_parzen_window_power.RData',
                normal='comparisons/out/NOrMAL.Rdata',
                local='powerAnalysis/output_detection_power_pm3.RData')
kOutputPath <- 'powerAnalysis/results_design_regressions.RData'
kGeneLength <- 3501

# Function definitions

LoadIntoNewEnv <- function(path) {
  env.data <- new.env(hash=TRUE)
  load(path, envir=env.data)
  return(env.data)
}


NormalizeToUnitInterval <- function(x) {
  if (is.character(x)) {
    x <- as.numeric(x)
  }
  (x - min(x)) / diff(range(x))
}


NormalizeRegressors <- function(dd) {
  if ('coverage' %in% names(dd))
    dd$coverage.adj <- NormalizeToUnitInterval(dd$coverage)
  
  if ('magnitude' %in% names(dd))
    dd$magnitude.adj <- NormalizeToUnitInterval(dd$magnitude)
  
  if ('eff.magnitude' %in% names(dd))
    dd$eff.magnitude.adj <- NormalizeToUnitInterval(dd$eff.magnitude)
  
  if ('offset' %in% names(dd))
    dd$offset.adj <- NormalizeToUnitInterval(dd$offset)
  
  if ('offset.primary' %in% names(dd))
    dd$offset.primary.adj <- NormalizeToUnitInterval(dd$offset.primary)
  
  return(dd)
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

# Create adjusted covariates for interpretability
# Normalizing to range 1
data.list <- lapply(data.list, NormalizeRegressors)

results.list$cluster$power.cluster <- NormalizeRegressors(
  results.list$cluster$power.cluster)
results.list$normal$power.cluster <- NormalizeRegressors(
  results.list$normal$power.cluster)
results.list$parzen$power.cluster <- NormalizeRegressors(
  results.list$parzen$power.cluster)
results.list$local$power.primary <- NormalizeRegressors(
  results.list$local$power.primary)
results.list$local$power.alternative <- NormalizeRegressors(
  results.list$local$power.alternative)

# Build linear models for errors
lm.cluster.error <- lm(
  mean.abs.error ~ coverage.adj * offset.primary.adj * eff.magnitude.adj,
  weights=w.hat,
  data=data.list$cluster)
lm.normal.error <- lm(
  mean.abs.error ~ coverage.adj * offset.primary.adj * eff.magnitude.adj,
  weights=w.hat,
  data=data.list$normal)
lm.parzen.error <- lm(
  mean.abs.error ~ coverage.adj * offset.primary.adj * eff.magnitude.adj,
  weights=w.hat,
  data=data.list$parzen)
lm.primary.error <- lm(
  mean.abs.error ~ coverage.adj * offset.primary.adj * eff.magnitude.adj,
  weights=w.hat,
  data=data.list$primary)
lm.alternative.error <- lm(
  mean.abs.error ~ coverage.adj * offset.adj * eff.magnitude.adj,
  weights=w.hat,
  data=data.list$alternative)

# Build generalized linear models for power
glm.cluster.power <- glm(
  p.detected.primary ~ coverage.adj  * offset.primary.adj * eff.magnitude.adj,
  family=binomial(link='logit'),
  weights=n,
  data=results.list$cluster$power.cluster)
glm.normal.power <- glm(
  p.detected.primary ~ coverage.adj  * offset.primary.adj * eff.magnitude.adj,
  family=binomial(link='logit'),
  weights=n,
  data=results.list$normal$power.cluster)
glm.parzen.power <- glm(
  p.detected.primary ~ coverage.adj  * offset.primary.adj * eff.magnitude.adj,
  family=binomial(link='logit'),
  weights=n,
  data=results.list$parzen$power.cluster)
glm.primary.power <- glm(
  p.detected.primary ~ coverage.adj  * offset.primary.adj * eff.magnitude.adj,
  family=binomial(link='logit'),
  weights=n,
  data=results.list$local$power.primary)
glm.alternative.power <- glm(
  p.detected.alternative ~ coverage.adj  * offset.adj  * eff.magnitude.adj,
  family=binomial(link='logit'),
  weights=n,
  data=results.list$local$power.alternative)

# Save regressions
save(lm.cluster.error, glm.cluster.power,
     lm.normal.error, glm.normal.power,
     lm.parzen.error, glm.parzen.power,
     lm.primary.error, glm.primary.power,
     lm.alternative.error, glm.alternative.power,
     file=kOutputPath)
