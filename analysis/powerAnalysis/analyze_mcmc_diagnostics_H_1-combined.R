# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(ascii)

source('powerAnalysis/lib.R')

# Constants
kSummariesPattern <- 'results/mcmc_summaries_H_1-combined_chrom%02d.txt'
kChromList <- 1:16
# kSummariesPattern <- 'results/mcmc_summaries_powerAnalysis_chrom%02d.txt'
# kChromList <- 1:10
kNEffField <- 7  # 7th column for effective sample size
kThetaPostmeanField <- 1  # 1st column for posterior mean of theta
kSummaryPath <- 'out/mcmc_diagnostics_H_1-combined.txt'
kPlotPath <- 'out/mcmc_diagnostics_H_1-combined.pdf'
kQuantiles <- c(0.0005, 0.005, 0.025, 0.05, 0.1,
                0.9, 0.95, 0.975, 0.995, 0.9995)

# Function definitions

ReadColumnViaPipe <- function(path, column, ...) {
  kPipePattern <- "awk '{print $%d}' %s" # Pattern for awk-based field extraction
  conn <- pipe(sprintf(kPipePattern, column, path))
  
  dat <- NULL
  tryCatch(dat <- scan(conn, ...), finally=close(conn))
  return(dat)
}


# Script ------------------------------------------------------------------

# Load data
n.eff.list <- lapply(kChromList, function(chrom)
  data.frame(chrom=chrom, n.eff=ReadColumnViaPipe(
    sprintf(kSummariesPattern, chrom), col=kNEffField, skip=1)))
n.eff.df <- Reduce(rbind, n.eff.list)
rm(n.eff.list)

theta.postmean.list <- lapply(kChromList, function(chrom)
  data.frame(chrom=chrom, theta.postmean=ReadColumnViaPipe(
    sprintf(kSummariesPattern, chrom), col=kThetaPostmeanField, skip=1)))
theta.postmean.df <- Reduce(rbind, theta.postmean.list)
rm(theta.postmean.list)


# Compute summary statistics
n.eff.summaries <- c(summary(n.eff.df$n.eff, SD=sd(n.eff.df$n.eff)),
                     quantile(n.eff.df$n.eff, kQuantiles))


# Plot construction -------------------------------------------------------


# Build nice histogram
n.eff.hist <- qplot(x=n.eff.df$n.eff, geom='histogram', binwidth=25,
                    xlab=expression(paste('Effective sample size for ',
                                          theta[k])),
                    ylab='Frequency')

# Build scatterplot of effective sample size vs. theta_postmean
n.eff.scatter <- ggplot(mapping=aes(x=theta.postmean.df$theta.postmean,
                                    y=n.eff.df$n.eff))
n.eff.scatter <- n.eff.scatter + stat_density2d(geom='tile', contour=FALSE,
                                                aes(fill=..density..))
n.eff.scatter <- n.eff.scatter + ylab(
  expression(paste('Effective sample size for ', theta[k])))
n.eff.scatter <- n.eff.scatter + xlab(
  expression(paste('Posterior mean of ', theta[k])))



# Output ------------------------------------------------------------------

# Write summary statistics to text file
print(ascii(data.frame(n.eff.summaries)), type='rest', file=kSummaryPath)

# Write plots to PDF
pdf(kPlotPath, width=5, height=3.1, onefile=TRUE)

theme.default <- theme_get()
theme_set(theme_bw())
print(n.eff.hist)
print(n.eff.scatter)
theme_set(theme.default)

dev.off()
