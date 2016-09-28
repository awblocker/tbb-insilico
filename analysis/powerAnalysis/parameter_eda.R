# Load libraries
library(plyr)
library(ggplot2)
library(stringr)

source('powerAnalysis/lib.R')

# Constants

kNChrom <- 16

# Parameters

# Data paths

# Index of ORFs
path.gene.index <- '../geneInfo/geneIndex.txt'

# Path for template
path.template <- 'templates/template_H_1-combined.txt'

# Path for regions
path.regions <- '../XuData/regionDefs/regions_H_1-combined_800.txt'

# Path for reads
path.reads <- '../XuData/H_1-combined_chrom.txt'

# Pattern for parameter information
pattern.param <- 'results/mcmc_params_H_1-combined_chrom%02d.txt'


# Load data
regions <- ReadRaggedArray(path.regions, what=integer(), sep=' ',
                           index=TRUE, value.var='region', row.var='chrom',
                           index.var='pos')
reads <- ReadRaggedArray(path.reads, what=double(), sep=',', index=TRUE,
                         value.var='y', row.var='chrom', index.var='pos')

# Load parameters
param.paths <- sprintf(pattern.param, 1:kNChrom)
param.list <- lapply(param.paths, read.table, header=TRUE)
for (i in 1:length(param.list))
    param.list[[i]]$chrom <- i
params <- ldply(param.list, identity)

# Fix variable names
params <- rename(params, c('region_id'='region'))
names(params) <- str_replace(names(params), fixed('_'), '.')
colnames(params) <- names(params)

# Compute additional parameter summaries
params$y.mean.postmed <- exp(params$mu.postmed + params$sigmasq.postmed / 2.)

# Add region information to regions data
reads$region <- regions$region

# Cleanup
rm(regions)
gc()

# Aggregate reads information by region
read.summaries <- ddply(reads, .(chrom, region), summarise,
                        coverage=mean(y), sd.y=sd(y))

# Compute bins of coverage by region
read.summaries$coverage.bin <- cut_number(read.summaries$coverage, 10)

# Merge parameter information with read summaries
read.summaries <- join(read.summaries, params, by=c('chrom', 'region'),
                       type='left')

# Summarise parameters by coverage
params.by.coverage <- ddply(read.summaries, .(coverage.bin), summarise,
                            mu.mean=mean(mu.postmean),
                            mu.med=median(mu.postmean),
                            sigma.mean=mean(sigma.postmean),
                            sigma.med=median(sigma.postmean))

