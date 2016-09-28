library(reshape2)
library(data.table)
library(ggplot2)
library(Biostrings)

kAlbert <- "../Albert2007/stats_yeast_H2AZ_2006_05_01.csv"
kInputs <- c(
  "titration-1"="stats_Sample_NucDG_7_CAGATC_chr01.bam_chr01.csv",
  "titration-2"="stats_Sample_NucDG_8_ACTTGA_chr01.bam_chr01.csv",
  "titration-3"="stats_Sample_NucDG_9_GATCAG_chr01.bam_chr01.csv",
  "titration-4"="stats_Sample_NucDG_10_TAGCTT_chr01.bam_chr01.csv")

SumAN <- function(x) {
  sum(as.numeric(x))
}

AddWeights <- function(.dt, alpha = 0) {
  .dt[, w := ((control_frequency + alpha) / SumAN(control_frequency + alpha)) /
      ((cut_frequency + alpha) / SumAN(cut_frequency + alpha))]
  return(.dt)
}

AddPct <- function(.dt, by = c()) {
  .dt[, pct_control := control_frequency / SumAN(control_frequency) * 100,
      by = by]
  .dt[, pct_cut := cut_frequency / SumAN(cut_frequency) * 100,
      by = by]
  return(.dt)
}

LoadAlbert <- function(path) {
  albert <- data.table(read.csv(kAlbert))
  albert <- AddWeights(albert)
  albert <- AddPct(albert)
  albert[, sample := "albert2007"]
  return(albert)
}

LoadStats <- function(path) {
  stats <- fread(path)
  stats[, cuts := paste(sort(c(dinucleotide1, dinucleotide2)), collapse=','),
        by=list(dinucleotide1, dinucleotide2)]
  setkey(stats, cuts)
  stats <- AddWeights(stats)
  return(stats)
}

GetMarginals <- function(stats) {
  marginal <- stats[, lapply(.SD, sum), keyby=cuts,
                    .SDcols=c("cut_frequency", "control_frequency")]
  AddWeights(marginals)
  return(marginals)
}

GetSingle <- function(stats) {
  long <- melt(stats,
               id.vars=c("cut_frequency", "control_frequency"),
               measure.vars=c("dinucleotide1", "dinucleotide2"),
               value.name="dinucleotide")
  single <- long[, lapply(.SD, sum), keyby = "dinucleotide",
                 .SDcols=c("cut_frequency", "control_frequency")]
  single <- AddPct(single)
  single <- AddWeights(single)
  return(single)
}

GetMargins <- function(stats) {
  margins <- lapply(c("1", "2"), function(end) {
                    key <- paste0("dinucleotide", end)
                    margin <- stats[, lapply(.SD, sum), keyby = key,
                                    .SDcols = c("cut_frequency",
                                                "control_frequency")]
                    margin[, end := end]
                    setnames(margin, key, "cut")
                    return(margin)
               })
  margins <- rbindlist(margins)
  margins <- AddPct(margins)
  return(margins)
}

Main <- function() {
  samples <- list()
  samples[["albert"]] <- LoadAlbert(kAlbert)
  for (f in names(kInputs)) {
    stats <- LoadStats(kInputs[[f]])
    s <- GetSingle(stats)
    s[, sample := f]
    samples[[f]] <- s
  }
  samples <- Reduce(function(x, y) rbind(x, y, use.names = TRUE), samples)
  setkey(samples, dinucleotide, sample)
  print(qplot(x = dinucleotide, y = pct_cut, fill = sample, data = samples,
        geom = "bar", stat = "identity", position = "dodge"))
  return(samples)
}

if (sys.nframe() == 0) {
  Main()
}
