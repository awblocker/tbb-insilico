# Load libraries
library(ggplot2)
library(plyr)
library(stringr)
library(reshape2)
library(xtable)

source('powerAnalysis/lib.R')

# Set constants
kInputPath <- 'powerAnalysis/results_design_anovas.RData'
kAnovaTablePath <- 'powerAnalysis/tables_design_anovas.tex'
kAodTablePath <- 'powerAnalysis/tables_design_aods.tex'


# Function definitions ----------------------------------------------------



BuildXtable <- function(name, title, ...) {
  xtable(get(name, envir=parent.frame()), caption=title, ...)
}

# Script ------------------------------------------------------------------

# Load data

load(kInputPath)

# Find lm and glm objects
names.anovas <- ls()[str_detect(ls(), perl('^anova\\.(.*)\\.error$'))]
titles.anovas <- str_replace(names.anovas, perl('^anova\\.(.*)\\.error$'),
                             '\\1')

names.aods <- ls()[str_detect(ls(), perl('^anova\\.(.*)\\.power$'))]
titles.aods <- str_replace(names.aods, perl('^anova\\.(.*)\\.power$'),
                           '\\1')

# Build xtable object for each lm and glm object loaded
xtables.anovas <- mapply(BuildXtable, names.anovas, titles.anovas, SIMPLIFY=FALSE)
xtables.aods <- mapply(BuildXtable, names.aods, titles.aods, SIMPLIFY=FALSE)

# Save xtables to appropriate tex files
for (ii in 1:length(xtables.anovas)) {
  print(xtables.anovas[[ii]], file=kAnovaTablePath, append=(ii > 1))
}

for (ii in 1:length(xtables.aods)) {
  print(xtables.aods[[ii]], file=kAodTablePath, append=(ii > 1))
}
