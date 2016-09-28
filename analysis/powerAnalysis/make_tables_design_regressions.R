# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(xtable)

source('powerAnalysis/lib.R')

# Set constants
kInputPath <- 'powerAnalysis/results_design_regressions.RData'
kLmTablePath <- 'powerAnalysis/tables_design_lm.tex'
kGlmTablePath <- 'powerAnalysis/tables_design_glm.tex'


# Function definitions


BuildXtable <- function(name, title, ...) {
  xtable(get(name, envir=parent.frame()), caption=title, ...)
}

# Script ------------------------------------------------------------------

# Load data

load(kInputPath)

# Find lm and glm objects
names.lms <- ls()[str_detect(ls(), perl('^lm\\.(?!names)'))]
titles.lms <- str_replace(names.lms, perl('lm\\.(.*)'), '\\1')

names.glms <- ls()[str_detect(ls(), perl('^glm\\.(?!names)'))]
titles.glms <- str_replace(names.glms, perl('glm\\.(.*)'), '\\1')

# Build xtable object for each lm and glm object loaded
xtables.lms <- mapply(BuildXtable, names.lms, titles.lms, SIMPLIFY=FALSE)
xtables.glms <- mapply(BuildXtable, names.glms, titles.glms, SIMPLIFY=FALSE)

# Save xtables to appropriate tex files
for (ii in 1:length(xtables.lms)) {
  print(xtables.lms[[ii]], file=kLmTablePath, append=(ii > 1))
}

for (ii in 1:length(xtables.glms)) {
  print(xtables.glms[[ii]], file=kGlmTablePath, append=(ii > 1))
}


