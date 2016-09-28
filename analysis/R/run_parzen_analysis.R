# Load libraries
library(plyr)
library(stringr)
library(reshape2)
library(yaml)
library(argparse)

source('R/lib.R')

# Set constants

kProg <- 'Rscript run_parzen_analysis.R'
kDescription <- 'Run Parzen window-based analysis.'
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
  '--tags', default='', type='character',
  help='Comma-separated list of tags to parse in CONFIG.')
kParser$add_argument(
  '--bw', default=20, type='double',
  help='Bandwidth (standard deviation) for window [default %(default)s]')
kParser$add_argument(
  '--wbw', default=3, type='double',
  help='Width (+/- standard deviations) for window [default %(default)s]')
kParser$add_argument(
  '--min-distance', default=147, type='double',
  help='Minimal distance between called peaks [default %(default)s]')


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
  
  # Load reads
  reads.list <- lapply(strsplit(readLines(cfg$data$chrom_path), ','),
                       as.integer)
  
  # Run Parzen window smoothing and peak calling
  parzen.list <- lapply(reads.list, CallPeaksParzenWindow,
                        bw=optargs$bw, w.bw=optargs$wbw,
                        min.distance=optargs$min_distance,
                        boundary=FALSE)
  
  # Extract smoothed results
  smoothed.df <- Reduce(rbind, mapply(function(parzen, chrom) {
    out <- data.frame(chrom=chrom,
                      pos=seq_along(parzen$y.smoothed),
                      smoothed=as.numeric(parzen$y.smoothed))
    return(out)
  }, parzen.list, optargs$chrom, SIMPLIFY=FALSE))
  
  # Extract peaks
  peaks.df <- Reduce(rbind, mapply(function(parzen, chrom) {
    out <- data.frame(chrom=chrom,
                      peak=parzen$peaks)
    return(out)
  }, parzen.list, optargs$chrom, SIMPLIFY=FALSE))
  
  # Save output
  cfg.name <- str_replace(basename(config.path), perl('\\.\\w*$'), '')
  
  smoothed.path <- file.path(
    optargs$outdir,sprintf('parzen_smoothed_%s_bw%d_dist%d.txt',
                           cfg.name, optargs$bw, optargs$min_distance))
  peaks.path <- file.path(
    optargs$outdir, sprintf('parzen_peaks_%s_bw%d_dist%d.txt',
                            cfg.name, optargs$bw, optargs$min_distance))
  
  write.table(smoothed.df, file=smoothed.path,
              row.names=FALSE, col.names=TRUE, quote=FALSE)
  write.table(peaks.df, file=peaks.path,
              row.names=FALSE, col.names=TRUE, quote=FALSE)
}

