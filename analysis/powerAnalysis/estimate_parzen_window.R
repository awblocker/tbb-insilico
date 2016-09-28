# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)

source('powerAnalysis/lib.R')

# Constants
kRegionPath <- 'powerAnalysis/data/simChrom_regions.txt'  # Path for regions
kReadPath <- 'powerAnalysis/data/simChrom_y.txt'  # Path for reads
kNullPath <- 'powerAnalysis/data/simChrom_null.txt'  # Path for null reads
kBw <- 20  # Bandwidth for Parzen window estimation
kWbw <- 3  # Width of window (+/-) as multiple of bw
kMinDistance <- 147 # Minimum distance between called clusters
kBoundary <- FALSE  # Do not include calls the boundaries of chromosomes
kPeakPath <- 'powerAnalysis/parzen_window_peaks.RData'
kNullPeakPath <- 'powerAnalysis/parzen_window_null_peaks.RData'


# Functions

RunParzenWindowAnalysis <- function(reads) {
  # Run peak calling on each chromosome
  peaks <- ddply(reads, .(chrom), function(df)
    data.frame(peaks=CallPeaksParzenWindow(y=df$y, bw=kBw, w.bw=kWbw,
                                           min.distance=kMinDistance,
                                           boundary=kBoundary)$peaks),
    .progress='text')
  return(peaks)
}


# Run Parzen window estimation for nonnull data
reads <- ReadRaggedArray(kReadPath, what=double(), sep=',', index=TRUE,
                         value.var='y', row.var='chrom', index.var='pos')
peaks <- RunParzenWindowAnalysis(reads)

# Save output
save(peaks, file=kPeakPath)

# Cleanup
rm(peaks, reads)
gc()

# Run Parzen window estimation for null data
null<- ReadRaggedArray(kNullPath, what=double(), sep=',', index=TRUE,
                       value.var='y', row.var='chrom', index.var='pos')
null.peaks <- RunParzenWindowAnalysis(null)

# Save output
save(null.peaks, file=kNullPeakPath)

