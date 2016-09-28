# Load libraries
library(plyr)
library(stringr)
library(reshape2)
library(yaml)
library(argparse)

source('R/lib.R')

# Set constants

kProg <- 'Rscript extract_detections.R'
kDescription <- 'Extract and combine detections.'
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
  '--pm', type='integer', default=1,
  help='+/- setting to extract [default %(default)s]')
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
  
  # Load detections
  cfg.name <- str_replace(basename(config.path), perl('\\.\\w*$'), '')
  
  detections.df <- mdply(
    data.frame(chrom=optargs$chrom),
    function(chrom, pm, pattern, ...)
      read.table(sprintf(pattern, chrom, pm), header=TRUE),
    pm=optargs$pm, pattern=cfg$mcmc_output$detections_pattern)
  
  # Save output
  
  out.path <- file.path(optargs$outdir,
                        sprintf('detections_%s_pm%d.txt', cfg.name, optargs$pm))
  write.table(detections.df, file=out.path, row.names=FALSE, col.names=TRUE,
              quote=FALSE)
}
