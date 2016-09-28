# Load libraries
library(argparse)
library(data.table)

library(stringr)
library(reshape2)
library(yaml)
library(ascii)

library(qvalue)
library(fdrtool)

options(stringsAsFactors = FALSE)

source('R/lib.R')

kProg <- 'Rscript analyze_local_control.R'
kDescription <- 'Run FDR analysis based on observed and control posterior summaries.'
kEpilogue <- 'Requires the argparse R package.\\n\\n'

kParser <- ArgumentParser(description = kDescription, epilog = kEpilogue,
                          prog = kProg)

kParser$add_argument(
  'config', metavar = 'CONFIG', type = 'character', nargs = 1,
  help = 'A YAML configuration file.')
kParser$add_argument(
  'outpath', metavar = 'OUTFILE', type = 'character', nargs = 1,
  help = 'Path for text output.')
kParser$add_argument(
  '--chrom', default = '', type = 'character',
  help = 'Comma-separated list of chromosomes to analyze. Defaults to all.')
kParser$add_argument(
  '--mean_regexp', type = 'character', default = 'mean_local_concentration.*',
  help = 'Regular expression (PCRE) for means [default %(default)s]')
kParser$add_argument(
  '--se_regexp', type = 'character', default = 'se_local_concentration.*',
  help = 'Regular expression (PCRE) for std errors [default %(default)s]')
kParser$add_argument(
  '--fdr', default = '0.2,0.1,0.05,0.01,0.005,0.001',
  help = 'Comma-separated list of FDRs to analyze [default %(default)s].')
kParser$add_argument(
  '--replicates', default = 100, type = 'integer',
  help = 'Number of replicates for permutation null [default %(default)s].')
kParser$add_argument(
  '--tags', default = '', type = 'character',
  help = 'Comma-separated list of tags to parse in CONFIG.')


# Function definitions

LoadViaList <- function(paths, ids, id.var = 'id', load.fun = fread, ...) {
  data.list <- lapply(paths, load.fun, ...)

  for (i in 1:length(data.list)) {
    data.list[[i]][[id.var]] <- ids[i]
  }

  data.df <- rbindlist(data.list)
  rm(data.list)
  gc()

  return(data.df)
}

LoadConfig <- function(optargs) {
  cfg <- yaml.load_file(optargs$config)
  if (str_length(optargs$tags[1]) > 0) {
    cfg <- ParseConfig(config = cfg, tags = optargs$tags)
  }
  if (length(optargs$chrom) == 0) {
    optargs$chrom <- 1:cfg$data$n_chrom
  }
  return(cfg)
}

TestStatistic <- function(experiment.mean, experiment.se,
                          control.mean, control.se) {
  # Using a two-sample z-statistic to get local testing along the genome
  # between an experimental condition and a control.
  return((experiment.mean - control.mean) /
         sqrt(experiment.se^2 + control.se^2))
}

EstimateTwoGroupModel <- function(stat) {
  # Parameterizing as alternative mean, log null SD,
  # log alternative SD / null SD, and p(null).
  negll <- function(theta) {
    f.null <- dt(stat / exp(theta[2]), theta[5]) / exp(theta[2])
    f.alt <- dt((stat - theta[1]) / exp(theta[2] + theta[3]), theta[5]) /
      exp(theta[2] + theta[3])
    -sum(log(theta[4] * f.null + (1 - theta[4]) * f.alt))
  }
  mle <- optim(c(mean(stat), log(sd(stat)), 0, 0.5, 30),
               negll,
               method = "L-BFGS-B",
               lower = c(-Inf, -Inf, -Inf, 0, 1),
               upper = c(Inf, Inf, Inf, 1, Inf))
  p.null <- mle$par[4]
  mean.alt <- mle$par[1]
  scale.null <- exp(mle$par[2])
  scale.alt <- exp(mle$par[2] + mle$par[3])
  df <- mle$par[5]
  return(function(x) p.null * dt(x / scale.null, df) / scale.null /
           (p.null * dt(x / scale.null, df) / scale.null +
              (1 - p.null) * dt((x - mean.alt) / scale.alt, df) / scale.alt))
}

GetFdrInformation <- function(obs.stat, null.cdf, fdrs, ...) {
  # Compute p-values of observed statistics
  obs.p.values <- 1 - null.cdf(obs.stat)
  stats <- list(
    Unadjusted = obs.p.values,
    qvalue = qvalue(obs.p.values)$qvalue,
    BH = p.adjust(obs.p.values, method = 'BH'),
    BY = p.adjust(obs.p.values, method = 'BY'))
  results <- list()
  # Iterate over FDRs of interest
  for (fdr in fdrs) {
    for (method in names(stats)) {
      info <- data.table(
        fdr = fdr, method = method,
        threshold.p.value = max(obs.p.values[stats[[method]] <= fdr]),
        threshold.stat = min(obs.stat[stats[[method]] <= fdr]))
      results[[length(results) + 1]] <- info
    }
  }
  return(rbindlist(results))
}

ToLong <- function(df) {
  df$position <- seq(nrow(df))
  melt(df, id.vars = c("chrom", "position"))
}

Main <- function(argv) {
  # Parse options
  optargs <- kParser$parse_args(argv)

  optargs$fdr <- as.numeric(strsplit(optargs$fdr, ',', fixed = TRUE)[[1]])
  optargs$chrom <- as.integer(strsplit(optargs$chrom, ',', fixed = TRUE)[[1]])

  # Load configuration file
  cfg <- LoadConfig(optargs)

  # Load posterior summaries
  obs <- LoadViaList(
    paths = sprintf(cfg$mcmc_output$summary_pattern, optargs$chrom),
    ids = optargs$chrom, id.var = 'chrom', load.fun = fread)

  control <- LoadViaList(
    paths = sprintf(cfg$mcmc_output$null_summary_pattern, optargs$chrom),
    ids = optargs$chrom, id.var = 'chrom', load.fun = fread)

  # Get names of input statistics and standard errors. Taking the intersection
  # here because we need each statistic to appear in both the observed and the
  # control summaries.
  names.stats <- data.table(
    mean.var = intersect(
      names(obs)[str_detect(names(obs), perl(optargs$mean_regexp))],
      names(control)[str_detect(names(obs), perl(optargs$mean_regexp))]),
    stderr.var = intersect(
      names(obs)[str_detect(names(obs), perl(optargs$se_regexp))],
      names(control)[str_detect(names(obs), perl(optargs$se_regexp))]),
    key = "mean.var")

  # Run analysis on each test statistic
  fdr.output <- list()
  for (mean.var in names.stats$mean.var) {
    stderr.var <- names.stats[mean.var]$stderr.var
    obs.stat <- TestStatistic(obs[[mean.var]], obs[[stderr.var]],
                              control[[mean.var]], control[[stderr.var]])
    null.cdf <- PermutationNull(control[[mean.var]], control[[stderr.var]],
                                num.replicates = optargs$replicates)
    fdr.output[[mean.var]] <- data.table(statistic = mean.var,
                                         GetFdrInformation(obs.stat = obs.stat,
                                                           null.cdf = null.cdf,
                                                           fdrs = optargs$fdr))
  }
  fdr.output <- rbindlist(fdr.output)

  # Save output
  print(ascii(fdr.output, include.rownames = FALSE, header = FALSE,
              format = 'g'),
        format = 'rest',
        file = optargs$outpath)
}

if (sys.nframe() == 0) {
#   argv <- commandArgs(TRUE)
  argv <- c("--chrom=1",
            "config/titration_control_2.yml",
            "out/fdr_titration_control_2.yml")
  Main(argv)
}
