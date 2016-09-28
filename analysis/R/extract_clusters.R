# Load libraries
library(plyr)
library(stringr)
library(reshape2)
library(yaml)
library(argparse)

source('R/lib.R')

# Set constants
kColumns <- c('center', 'occupancy', 'occupancy_se')
kProg <- 'Rscript extract_clusters.R'
kDescription <- 'Extract and combine clusters.'
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
  
  # Load clusters
  cfg.name <- str_replace(basename(config.path), perl('\\.\\w*$'), '')
  
  clusters.df <- mdply(
    data.frame(chrom=optargs$chrom),
    function(chrom, pattern, ...)
      read.table(sprintf(pattern, chrom), header=TRUE)[, kColumns],
    pattern=cfg$mcmc_output$cluster_pattern)
  
  # Save output
  
  out.path <- file.path(optargs$outdir,
                        sprintf('clusters_%s.txt', cfg.name))
  write.table(clusters.df, file=out.path, row.names=FALSE, col.names=TRUE,
              quote=FALSE)
}
