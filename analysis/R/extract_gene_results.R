# Load libraries
library(plyr)
library(reshape2)
library(stringr)
library(yaml)

source('R/lib.R')

# Set constants

# kGenes <- c("PHO11")
kGenes <- c("PHO11", "PHO12", "PHO13", "PHO2", "PHO23", "PHO3", "PHO4", "PHO5", 
            "PHO8", "PHO80", "PHO81", "PHO84", "PHO85", "PHO86", "PHO87",
            "PHO88", "PHO89", "PHO90", "PHO91", "MSN2", "MSN4")
kOutputDir <- "genes"

kConfigPath <- 'config/titration_1.yml'
kGeneIndexPath <- '../geneInfo/geneIndex.txt'
kDetectionPm <- 1
kPromoter <- 1000
kGeneSubset <- seq(-500, 500)
kBw <- 20  # Bandwidth for Parzen window estimation
kWbw <- 3  # Width of window (+/-) as multiple of bw
kMinDistance <- 147 # Minimum distance between called clusters
kBoundary <- FALSE

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


GetGenePositions <- function(common, gene.index, promoter=1000, subset=NULL){ 
  # Build vector of positions in gene
  gene <- gene.index[gene.index$common==common,]
  chrom <- gene$chrom
  if (gene$start < gene$stop) {
    inGene <- seq(gene$start-promoter, gene$stop)
  } else {
    inGene <- seq(gene$stop, gene$start+promoter)
  }
  
  # Setup vector of offsets from TSS
  relPos <- inGene - gene$start
  if (gene$stop < gene$start)
    relPos <- -relPos
  
  gene.df <- data.frame(chrom=chrom, pos=inGene, rel.pos=relPos)
  
  if (!is.null(subset))
    gene.df <- gene.df[gene.df$rel.pos %in% subset, ]
  
  return(gene.df)
}


GetReads <- function(cfg, gene.df, reads.list=NULL) {
  if (is.null(reads.list))
    reads.list <- lapply(strsplit(readLines(cfg$data$chrom_path), ','),
                         as.integer)
  
  reads.df <- gene.df
  reads.df$y <- reads.list[[gene.df$chrom[1]]][gene.df$pos]
  reads.df <- reads.df[order(reads.df$rel.pos), ]
  
  return(reads.df)
}


GetDetections <- function(cfg, gene.df, detection.pm) {
  # Load detections
  detections <- read.table(sprintf(cfg$mcmc_output$detections_pattern,
                                   gene.df$chrom[1], detection.pm), header=TRUE)
  
  # Subset detections to gene
  detections <- detections[detections$pos %in% gene.df$pos, ]
  
  # Add relative positions and sort by them
  detections$rel.pos <- gene.df$rel.pos[match(detections$pos, gene.df$pos)]
  detections <- detections[order(detections$rel.pos), ]
  detections$chrom <- gene.df$chrom[1]
  
  return(detections)
}


GetCoef <- function(cfg, gene.df) {
  summaries.path <- sprintf(cfg$mcmc_output$summary_pattern, gene.df$chrom[1])
  all.df <- LoadColumnsFromRegexp(path=summaries.path, regexp.var='^(se_)?b$')
  
  b.df <- gene.df
  b.df <- data.frame(b.df, all.df[gene.df$pos,])
  b.df <- b.df[order(b.df$rel.pos), ]
  
  return(b.df)
}


GetClusters <- function(cfg, gene.df) {
  clusters.path <- sprintf(cfg$mcmc_output$cluster_pattern, gene.df$chrom[1])
  all.df <- read.table(file=clusters.path, header=TRUE)
  
  clusters.df <- all.df[all.df$center %in% gene.df$pos, ]
  clusters.df$rel.pos <- gene.df$rel.pos[match(clusters.df$center, gene.df$pos)]
  clusters.df$chrom <- gene.df$chrom[1]
  clusters.df <- clusters.df[order(clusters.df$rel.pos), ]
  
  return(clusters.df)
}


GetParzenResults <- function(parzen.list, gene.df) {
  parzen.df <- gene.df
  parzen.df$smoothed <- parzen.list[[parzen.df$chrom[1]]]$y.smoothed[
    parzen.df$pos]
  parzen.df$is.peak <- as.integer(
    parzen.df$pos %in% parzen.list[[parzen.df$chrom[1]]]$peaks)
  
  return(parzen.df)
}


# Script

# Load configuration
cfg <- yaml.load_file(kConfigPath)

# Load gene index
gene.index <- read.table(kGeneIndexPath, header=TRUE)

# Load reads
reads.list <- lapply(strsplit(readLines(cfg$data$chrom_path), ','), as.integer)

# Run Parzen window smoothing and peak calling
parzen.list <- lapply(reads.list, CallPeaksParzenWindow,
                      bw=kBw, w.bw=kWbw, min.distance=kMinDistance,
                      boundary=kBoundary)

for (gene in kGenes) {
  # Extract data
  gene.df <- GetGenePositions(common=gene, gene.index=gene.index,
                              promoter=kPromoter, subset=kGeneSubset)
  reads <- GetReads(cfg=cfg, gene.df=gene.df, reads.list=reads.list)
  detections <- GetDetections(cfg=cfg, gene.df=gene.df,
                              detection.pm=kDetectionPm)
  b <- GetCoef(cfg=cfg, gene.df=gene.df)
  clusters <- GetClusters(cfg=cfg, gene.df=gene.df)
  parzen <- GetParzenResults(parzen.list=parzen.list, gene.df=gene.df)
  
  # Save data to zip archive
  files.vec <- paste(tempdir(),
                     c("reads.txt", "detections.txt", "b.txt", "clusters.txt",
                       "parzen.txt"),
                     sep=.Platform$file.sep)
  names(files.vec) <- c("reads", "detections", "b", "clusters", "parzen")
  
  for (nme in names(files.vec)) {
    write.table(get(nme), file=files.vec[nme], col.names=TRUE, row.names=FALSE)
  }
  
  zip.path <- paste(kOutputDir, sprintf("%s.zip", gene), sep=.Platform$file.sep)
  zip(zipfile=zip.path, files=files.vec, flags="-r9Xj")
}
