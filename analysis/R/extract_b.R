# Load libraries
library(plyr)
library(stringr)
library(reshape2)
library(yaml)
library(argparse)

source('R/lib.R')

# Set constants

kProg <- 'Rscript extract_b.R'
kDescription <- 'Extract posterior means and std. errors of beta.'
kEpilogue <- 'Requires the argparse R package.\\n\\n'

kParser <- ArgumentParser(description=kDescription, epilog=kEpilogue,
                          prog=kProg)

kParser$add_argument(
  'config', metavar='CONFIG', type='character', nargs='+',
  help='A YAML configuration file.')
kParser$add_argument(
  '--outdir', metavar='OUTDIR', type='character', default='.',
  help='Path to the output directory.')
kParser$add_argument(
  '--chrom', default='', type='character',
  help='Comma-separated list of chromosomes to analyze. Defaults to all.')
kParser$add_argument(
  '--regexp', type='character', default='^(se_)?b$',
  help='Regular expression (PCRE) for statistics [default %(default)s]')
kParser$add_argument(
  '--tags', default='', type='character',
  help='Comma-separated list of tags to parse in CONFIG.')


# Function definitions

LoadColumnsFromRegexp <- function(path, regexp.var, sep=perl('\\s')) {
  conn <- file(path)
  var.names <- readLines(conn, n=1)
  close(conn)
  
  var.names <- str_split(var.names, sep)[[1]]
  columns <- which(str_detect(var.names, regexp.var))
  var.names <- var.names[columns]
  
  column.list <- lapply(columns, function(col)
    ReadColumnViaPipe(path, col, skip=1))
  names(column.list) <- var.names
  return(as.data.frame(column.list))
}

LoadViaList <- function(paths, ids, id.var='id', load.fun=read.table,
                        index.var=NULL, ...) {
  data.list <- lapply(paths, load.fun, ...)
  
  for (i in 1:length(data.list)) {
    data.list[[i]][[id.var]] <- ids[i]
    if (!is.null(index.var))
      data.list[[i]][index.var] <- seq(nrow(data.list[[i]]))
  }
  
  data.df <- Reduce(rbind, data.list)
  rm(data.list)
  gc()
  
  return(data.df)
}


# Script

# Parse options
argv <- commandArgs(TRUE)
optargs <- kParser$parse_args(argv)

optargs$chrom <- as.integer(strsplit(optargs$chrom, ',', fixed=TRUE)[[1]])

for (config.path in optargs$config) {
  # Load configuration file
  cfg <- yaml.load_file(config.path)
  if (str_length(optargs$tags[1]) > 0)
    cfg <- ParseConfig(config=cfg, tags=optargs$tags)
  
  if (length(optargs$chrom) == 0)
    optargs$chrom <- 1:cfg$data$n_chrom
  
  # Load test statistics
  obs.df <- LoadViaList(
    paths=sprintf(cfg$mcmc_output$summary_pattern, optargs$chrom),
    ids=optargs$chrom, id.var='chrom', index.var='pos',
    load.fun=LoadColumnsFromRegexp, regexp.var=perl(optargs$regexp))
  
  # Save output
  cfg.name <- str_replace(basename(config.path), perl('\\.\\w*$'), '')
  out.path <- file.path(optargs$outdir, sprintf('mcmc_coef_%s.txt', cfg.name))
  write.table(obs.df, file=out.path, row.names=FALSE, col.names=TRUE,
              quote=FALSE)
}
