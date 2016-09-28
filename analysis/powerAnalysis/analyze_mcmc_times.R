# Load libraries
library(plyr)
library(reshape2)
library(ggplot2)
library(ascii)

# Set constants
kLogPaths <- list(sim='dump/mcmcPowerAnalysis.err',
                  H1='dump/mcmc_H_1-combined.err',
                  H2='dump/mcmc_H_2-combined.err')
kExtractionCommand <- "perl -nle 'if (/%s/) {print %s}' %s"
kLogRegexp <- '(\\d*):\\s*Iteration time: ([0-9.]*)'
kDataPattern <- '$1, ",", $2'
kSep <- ','
kDataType <- list(integer(0), numeric(0))
kColNames <- c('iteration', 'time')
kNIterations <- 2000


# Function definitions

ReadLog <- function(path, extraction.command, log.regexp, data.pattern,
                    data.type, col.names=NULL, ...) {
  # Open pipe to extract data
  conn <- pipe(sprintf(kExtractionCommand, kLogRegexp, kDataPattern, path))
  tryCatch(dataset <- scan(conn, data.type, ...), finally=close(conn))
  
  dataset <- as.data.frame(dataset)
  if (!is.null(col.names))
    names(dataset) <- colnames(dataset) <- col.names
    
  return(dataset)
}


# Script ------------------------------------------------------------------

# Load logs
log.list <- lapply(kLogPaths, ReadLog, extraction.command=kExtractionCommand,
                   log.regexp=kLogRegexp, data.pattern=kDataPattern,
                   data.type=kDataType, sep=kSep,
                   col.names=kColNames)

# Compress high-phosphate runs
log.HP <- log.list$H1
log.HP$time <- log.HP$time / 2 + log.list$H2$time / 2

# Compute summaries of runtimes
time.per.iteration.HP <- sapply(1:(nrow(log.HP) / kNIterations), function(i)
  mean(log.HP$time[(1 + (i-1)*kNIterations):(i*kNIterations)]))
time.per.chromosome.HP <- time.per.iteration.HP * kNIterations / 60

time.per.iteration.sim <- sapply(
  1:(nrow(log.list$sim) / kNIterations), function(i)
  mean(log.list$sim$time[(1 + (i-1)*kNIterations):(i*kNIterations)]))
time.per.chromosome.sim <- time.per.iteration.sim * kNIterations / 60
