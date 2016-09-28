# Load libraries
library(plyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(yaml)

source('powerAnalysis/lib.R')

# Set parameters

kConfigPath <- 'powerAnalysis/powerAnalysis.yml'
kSep <- ' '
kNRows <- -1
kVariableRegexp <- 'p_local_concentration.*'
kQuantiles <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.995, .999)
kThresholds <- seq(0.5, 0.95, 0.05)
kDataPath <- 'powerAnalysis/obs_mcmc_summaries.Rdata'
kOutputPath <- 'powerAnalysis/obs_distribution_summaries.Rdata'


# Script

# Load and parse configuration
config <- yaml.load_file(kConfigPath)
config <- ParseConfig(config)

# Load null results
obs.result.pattern <- config$mcmc_output$summary_pattern
obs.result.paths <- sprintf(obs.result.pattern,
                            1:as.integer(config$data$n_chrom))

obs.result.header <- scan(obs.result.paths[1], character(0), nlines=1)
obs.result.columns <- which(str_detect(obs.result.header, kVariableRegexp))
obs.result.colnames <- obs.result.header[obs.result.columns]
names(obs.result.columns) <- obs.result.colnames

obs.results.list <- lapply(obs.result.paths, function(path)
  data.frame(sapply(obs.result.columns, function(col)
    ReadColumnViaPipe(path, col, skip=1))))
for (i in 1:length(obs.results.list)) {
  obs.results.list[[i]]$replicate <- i
}

# Merge null results
obs.df <- Reduce(rbind, obs.results.list)

# Clean up
rm(obs.results.list)
gc()

# Find variables of interest for summaries
variables.of.interest <- names(obs.df)[str_detect(names(obs.df),
                                                  kVariableRegexp)]

# Compute tail probabilities and quantiles by replicate
tail.probs <- ddply(obs.df, .(replicate), ExtractTailProbabilities,
                    variables=variables.of.interest, thresholds=kThresholds)
quantiles <- ddply(obs.df, .(replicate), ExtractQuantiles,
                   variables=variables.of.interest, p=kQuantiles)

# Compute aggregate summaries of quantiles and tail probabilities
tail.prob.means <- dcast(melt(tail.probs, id.vars='threshold'),
                         threshold~variable, mean)
quantile.means <- dcast(melt(quantiles, id.vars='quantile'),
                        quantile~variable, mean)

# Save output
save(obs.df, file=kDataPath)
save(tail.probs, quantiles, tail.prob.means, quantile.means, file=kOutputPath)

