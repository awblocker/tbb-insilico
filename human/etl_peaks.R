library(rtracklayer)
library(ChIPpeakAnno)
library(AnnotationHub)

kChromosomes <- "chr20"
kPeaksUrl <- paste0("http://eqtl.uchicago.edu/nucleosomes/positioning_scores/",
                    "peaks.min_peak_score_0.6.thresh_0.5.txt.gz")
kArchivePath <- "output/peaks.Rdata"

ImportPeaks <- function(url) {
  tmp <- file.path(tempdir(), basename(url))
  download.file(url, tmp, method = "wget")
  peaks <- import.wig(tmp)
  peaks$score <- peaks$NA.
  peaks$NA. <- NULL
  return(peaks)
}

FilterPeaks <- function(peaks, chromosomes) {
  peaks[seqnames(peaks) %in% chromosomes, ]
}

AnnotatePeaks <- function(peaks, hub) {
  # RefSeq hg18 genes
  genes <- hub[["AH5155"]]
  annotated <- annotatePeakInBatch(RangedData(peaks),
                                   AnnotationData = genes,
                                   output = "both")
  annotated[["gene"]] <- genes$name[as.integer(annotated[["feature"]])]
  return(annotated)
}

Main <- function() {
  hub <- AnnotationHub()
  peaks <- ImportPeaks(kPeaksUrl)
  peaks <- FilterPeaks(peaks, kChromosomes)
  annotated <- AnnotatePeaks(peaks, hub)
  overlapping <- with(annotated,
                      annotated[insideFeature %in% c(
                        "inside", "overlapEnd", "overlapStart"), ])
  save(annotated, file = kArchivePath)
}

if (sys.nframe() == 0) {
  Main()
}
