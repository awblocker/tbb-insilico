library(Rcpp)

sourceCpp('powerAnalysis/lib.cpp')

# Functions

ReadColumnViaPipe <- function(path, column, ...) {
  kPipePattern <- "awk '{print $%d}' %s" # Pattern for awk-based field extraction
  conn <- pipe(sprintf(kPipePattern, column, path))
  
  dat <- NULL
  tryCatch(dat <- scan(conn, ...), finally=close(conn))
  return(dat)
}

#' Reads data from a ragged array, stored as delimited text
#'
#' @param path Path to delimited text file
#' @param what Atomic type
#' @param sep Separator for delimited text file
#' @param id Arbitrary identifier to add as extra column. Default is NULL.
#' @param index If TRUE, includes index along ragged rows. Default is FALSE.
#' @param value.var Name of variable to return containing values. Default is 
#'   'value'.
#' @param row.var Name of variable to return as row index. Default is 'row'.
#' @param index.var Variable name for index. Default is 'index'.
#' @param id.var Variable name for id. Default is 'id'.
#'
#' @return data.frame containing the row indices, values, and (optionally) index
#'   and id as columns
ReadRaggedArray <- function(path, what=double(), sep=' ', id=NULL, index=FALSE,
                            value.var='value', row.var='row', index.var='index',
                            id.var='id') {
    
    # Open connection to text file
    f <- file(path, 'rt')
    
    # Load text data
    rows <- list()
    while (length(line <- readLines(f, n=1)) > 0) {
        rows[[length(rows) + 1]] <- scan(text=line, what=what, sep=sep)
    }
    
    # Cleanup
    close(f)
    gc()

    # Arrange data into long-form data.frame

    # Allocate data.frame
    n.values <- sum(sapply(rows, length))
    data.long <- data.frame(row=vector('integer', n.values),
                            value=vector(storage.mode(what), n.values))
    if (index)
        data.long$index <- vector('integer', n.values)

    # Fill data.frame
    n.filled <- 0
    for (i in 1:length(rows)) {
        row.length <- length(rows[[i]])
        slice <- (n.filled + 1):(n.filled + row.length)
        data.long$row[slice] <- i
        data.long$value[slice] <- rows[[i]]
        if (index)
            data.long$index[slice] <- 1:row.length
        n.filled <- n.filled + row.length
    }
    
    # Add id, if requested
    if (!is.null(id))
        data.long$id <- as.factor(id)
    
    # Rename variables
    var.names <- c(row.var, value.var)

    if (index)
        var.names <- c(var.names, index.var)

    if (!is.null(id))
        var.names <- c(var.names, id.var)

    names(data.long) <- colnames(data.long) <- var.names

    gc()
    return(data.long)
}


MatchPositions <- function(x, y) {
    # Find nearest position in y for each position in x
    candidate.indices <- findInterval(x, y, all.inside=FALSE)
    n.x <- length(x)
    n.y <- length(y)

    candidate.positions <- matrix(NA, n.x, 2)
    candidate.positions[, 1] <- y[pmax(candidate.indices, 1)]
    candidate.positions[, 2] <- y[pmin(candidate.indices + 1, n.y)]

    candidate.distances <- abs(x - candidate.positions)

    matched.positions <- candidate.positions[
        cbind(1:nrow(candidate.positions),
              1 + (candidate.distances[, 1] > candidate.distances[, 2]))]

    return(matched.positions)
}


#' Find local maxima in a sequence
#' 
#' Simple definition of local maxima for a discrete sequence (low-high-low).
#' 
#' @param x Vector containing sequence to scan for local maxima
#' @param boundary Logical; include boundary cases? Default is FALSE.
#' @return vector containing integer indices of local maxima from sequence
FindLocalMaxima <- function(x, boundary=FALSE) {
  # Simple low-high-low definition
  left.gt.right <- c((head(x, -1) >= tail(x, -1)), boundary)
  right.gt.left <- c(boundary, (head(x, -1) <= tail(x, -1)))
  
  return(which(left.gt.right & right.gt.left))
}

#' Standard Parzen window estimation for nucleosome cluster positions
#' 
#' Runs field-standard Parzen window estimation for nucleosome cluster
#' postitions as described in Tirosh (2012).
#' 
#' @param y Numeric vector of read counts per base pair
#' @param bw Bandwidth for Gaussian kernel (standard deviation).
#' @param w.bw Width of kernel in multiples of bw (+/-).
#' @param min.distance Minimum distance between called peaks.
#' @param boundary Logical; include boundary cases? Default is FALSE.
#' @return List with three components:
#'  \itemize{
#'    \item y.smoothed, the smoothed read counts
#'    \item candidates, the local maxima of the smoothed read counts
#'    \item peaks, the filtered local peaks
#'  }
CallPeaksParzenWindow <- function(y, bw=25, w.bw=2, min.distance=120,
                                  boundary=FALSE) {
  # Build smoothing window
  w <- ceiling(bw * w.bw)
  window.gaussian <- dnorm(seq(-w, w), sd=bw)
  window.gaussian <- window.gaussian / sum(window.gaussian)
  
  # Smooth read counts
  y.smoothed <- filter(y, window.gaussian, method='convolution')
  y.smoothed[1:w] <- y.smoothed[w + 1]
  y.smoothed[(length(y) - w + 1):length(y)] <- y.smoothed[(length(y) - w)]
  
  # Find candidate peaks (local maxima)
  candidates <- FindLocalMaxima(y.smoothed, boundary=boundary)
  values.candidates <- y.smoothed[candidates]
  
  # Sort candidate peaks by values
  candidates <- candidates[order(values.candidates, decreasing=TRUE)]
  values.candidates <- sort(values.candidates)
  
  # Find peaks via greedy search
  peaks <- integer(length(candidates))
  n.peaks <- GreedySearch(candidates, min.distance, peaks)
  
  # Sort peaks for clarity
  peaks <- head(peaks, n.peaks)
  peaks <- sort(peaks)
  
  # Setup return value
  retval <- list(y.smoothed=y.smoothed, candidates=candidates, peaks=peaks)
  return(retval)
}

ParseConfig <- function(config, tags='id', value.tag=NULL) {
  # Parses simple Python-style string substitutions in config list. Runs
  # recursively to parse entire nested list structure.
  #
  # Args:
  #   config: A list generated from YAML config file
  #   tags: A character vector of substitutions to parse
  #   value.tag: Used in recursions. Leave as NULL.
  #
  # Returns:
  #   A list with the same structure as config with simple string
  #   substitutions parsed
  #
  # Run sequentially, once per tag
  for (tag in tags) {
    if (is.null(value.tag))
      value.tag <- config[[tag]]
    
    config <- lapply(config, function(entry) {
      if (is.list(entry))
        ParseConfig(entry, tags=tags, value.tag=value.tag)
      else
        str_replace(entry, fixed(str_c('{', tag, '}')),
                    value.tag)
    })
  }
  return(config)
}

ExtractQuantiles <- function(df, variables, p) {
  # Extract p quantiles of variables from data.frame.
  #
  # Args:
  #   df: data.frame containing variables as columns
  #   variables: variables to summarize
  #   p: quantiles to extract for each variable
  #
  # Returns:
  #   A data.frame with the requested quantiles.
  
  # Initialize for output
  quantiles <- data.frame(replicate(length(variables) + 1,
                                    numeric(length(p))))
  names(quantiles) <- c('quantile', variables)
  colnames(quantiles) <- names(quantiles)
  
  quantiles$quantile <- p
  
  # Iterate over variables
  for (v in variables) {
    quantiles[[v]] <- quantile(df[[v]], p)
  }
  return(quantiles)
}

ExtractTailProbabilities <- function(df, variables, thresholds) {
  # Extract P(x >= thresholds) for variables from data.frame.
  #
  # Args:
  #   df: data.frame containing variables as columns
  #   variables: variables to summarize
  #   p: quantiles to extract for each variable
  #
  # Returns:
  #   A data.frame with the requested empirical tail probabiltiies.
  
  # Initialize for output
  tail.probs <- data.frame(replicate(length(variables) + 1,
                                     numeric(length(thresholds))))
  names(tail.probs) <- c('threshold', variables)
  colnames(tail.probs) <- names(tail.probs)
  
  tail.probs$threshold <- thresholds
  
  # Iterate over variables
  for (v in variables) {
    tail.probs[[v]] <- sapply(thresholds, function(tr) mean(df[[v]] >= tr))
  }
  return(tail.probs)
}

CoefSe.glm <- function(obj) {
  beta.hat <- coef(obj)
  se.beta.hat <- sqrt(diag(summary(obj)$cov.scaled))
  return(data.frame(coef=beta.hat, se=se.beta.hat))
}
