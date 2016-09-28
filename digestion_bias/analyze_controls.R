# Validation for digestion-biased controls.
library(reshape2)
library(data.table)
library(ggplot2)
library(stringr)
library(knitr)
library(rmarkdown)

kDirPattern <- "control_Sample_NucDG*"
kReportFile <- "validate_controls.Rmd"

LoadExperiment <- function(dir) {
  cuts <- fread(file.path(dir, "cuts.csv"))
  experiment <- str_replace(dir, ignore.case("control_"), "")
  cuts[, experiment := experiment]
  lengths <- data.table(read.csv(file.path(dir, "lengths.csv"),
                                 col.names = c("length", "n")))
  lengths[, experiment := experiment]
  return(list(cuts = cuts, lengths = lengths))
}

LoadExperiments <- function(dirs) {
  experiments <- lapply(dirs, LoadExperiment)
  return(Reduce(function(a, b) list(cuts = rbind(a$cuts, b$cuts),
                                    lengths = rbind(a$lengths, b$lengths)),
                experiments))
}

CheckSymmetry <- function(cuts) {
  margins <- merge(cuts[, lapply(.SD, sum),
                        keyby = list(experiment, dinucleotide = fwd),
                        .SDcols = "num_reads"],
                   cuts[, lapply(.SD, sum),
                        keyby = list(experiment, dinucleotide = rev),
                        .SDcols = "num_reads"],
                   suffixes = c("_fwd", "_rev"))
  margins[, `:=`(
    lcb = qbeta(0.025, num_reads_fwd + 1, num_reads_rev + 1),
    ucb = qbeta(0.975, num_reads_fwd + 1, num_reads_rev + 1))]
  margins[, `:=`(
    lcb = lcb / (1 - lcb),
    ucb = ucb / (1 - ucb))]
  return(margins)
}

WriteReport <- function(experiments, symmetry) {
  rmarkdown::render(kReportFile)
}

RunningWithinRMarkdown <- function() {
  calls <- sys.calls()
  any(sapply(calls, function(v)
    any(str_detect(as.character(v), "rmarkdown"))))
}

Main <- function() {
  dirs <- list.files(".", pattern = kDirPattern)
  experiments <- LoadExperiments(dirs)
  symmetry <- CheckSymmetry(experiments$cuts)
  if (!RunningWithinRMarkdown()) {
    WriteReport(experiments = experiments, symmetry = symmetry)
  }
  return(list(experiments = experiments,
              symmetry = symmetry))
}

if (sys.nframe() == 0) {
  Main()
}
