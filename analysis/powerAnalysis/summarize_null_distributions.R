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
kDataPath <- 'powerAnalysis/null_mcmc_summaries.Rdata'
kOutputPath <- 'powerAnalysis/null_distribution_summaries.Rdata'

# Script

# Load and parse configuration
config <- yaml.load_file(kConfigPath)
config <- ParseConfig(config)

# Load null results
null.result.pattern <- config$mcmc_output$null_summary_pattern
null.result.paths <- sprintf(null.result.pattern,
                             1:as.integer(config$data$n_chrom))

null.result.header <- scan(null.result.paths[1], character(0), nlines=1)
null.result.columns <- which(str_detect(null.result.header, kVariableRegexp))
null.result.colnames <- null.result.header[null.result.columns]
names(null.result.columns) <- null.result.colnames

null.results.list <- lapply(null.result.paths, function(path)
  data.frame(sapply(null.result.columns, function(col)
    ReadColumnViaPipe(path, col, skip=1))))
for (i in 1:length(null.results.list)) {
    null.results.list[[i]]$replicate <- i
}

# Merge null results
null.df <- Reduce(rbind, null.results.list)

# Clean up
rm(null.results.list)
gc()

# Find variables of interest for summaries
variables.of.interest <- names(null.df)[str_detect(names(null.df),
                                                   kVariableRegexp)]

# Compute tail probabilities and quantiles by replicate
tail.probs <- ddply(null.df, .(replicate), ExtractTailProbabilities,
                    variables=variables.of.interest, thresholds=kThresholds)
quantiles <- ddply(null.df, .(replicate), ExtractQuantiles,
                   variables=variables.of.interest, p=kQuantiles)

# Compute aggregate summaries of quantiles and tail probabilities
tail.prob.means <- dcast(melt(tail.probs, id.vars='threshold'),
                         threshold~variable, mean)
quantile.means <- dcast(melt(quantiles, id.vars='quantile'),
                        quantile~variable, mean)

# Save output
save(null.df, file=kDataPath)
save(tail.probs, quantiles, tail.prob.means, quantile.means, file=kOutputPath)

