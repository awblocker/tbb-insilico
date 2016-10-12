library(data.table)
library(rtracklayer)

kMidpointsFile <- "/Users/awblocker/data/geo/chr20_combined_126_to_184.bed.gz"
kColumnNames <- c("chrom",
                  "start",
                  "end",
                  "score",
                  "geneChrom",
                  "geneStart",
                  "geneEnd",
                  "geneName")
kCountsFile <- "output/counts_chr20.txt"
kGenesFile <- "output/genes_chr20.txt"

ReadMidpoints <- function(filename) {
  midpoints <- fread(paste0("gzcat ", filename),
                     sep = "\t", header = FALSE,
                     select = c(1, 2, 3, 5, 6, 7, 8, 9))
  setnames(midpoints, kColumnNames)
  return(midpoints)
}

ConvertToCounts <- function(midpoints) {
  midpoints[, relativePosition := start - geneStart + 1]
  setkey(midpoints, geneName)
  geneNames <- unique(midpoints[, geneName])
  counts <- list()
  for (geneName in geneNames) {
    geneMidpoints <- midpoints[geneName]
    geneCounts <- integer(geneMidpoints[, geneEnd - geneStart + 1][1])
    geneCounts[geneMidpoints[, relativePosition]] <- geneMidpoints[, score]
    counts[[geneName]] <- geneCounts
  }
  return(counts)
}

WriteCounts <- function(counts, filename) {
  outfile <- file(filename, "wt")
  on.exit(close(outfile))
  for (gene in names(counts)) {
    writeLines(paste(counts[[gene]], collapse = ","), con = outfile)
  }
}

Main <- function() {
  midpoints <- ReadMidpoints(kMidpointsFile)
  counts <- ConvertToCounts(midpoints)
  genes <- unique(midpoints[, .(geneChrom, geneStart, geneEnd),
                            keyby = geneName])
  WriteCounts(counts, kCountsFile)
  write.csv(genes, kGenesFile, row.names = FALSE)
}

if (sys.nframe() == 0) {
  Main()
}
