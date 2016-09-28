# Load libraries
library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)

# Set parameters
templatePath    <- "templates/template_H_1-combined.txt"
lengthDistPath  <- "../XuData/H_1-combined_lengthDist.txt"
plotPathBase    <- "plots/plotTemplateVsApprox_H_1-combined"
convertCmdFmt   <- "pdftops -eps %s.pdf"
set.seed(19871013)

# Load template
template    <- scan(templatePath, 0.0)
w           <- length(template) %/% 2
offset      <- seq(-w,w)

# Get discrete Gaussian approximation
lengthDist  <- read.table(lengthDistPath, header=FALSE,
                            col.names=c("l","n"))
#
p           <- lengthDist$n/sum(lengthDist$n)
mu          <- sum(p*lengthDist$l)
sigma       <- sqrt( sum(p*(lengthDist$l-mu)^2) )
#
approx      <- exp(-(offset/sigma*2)^2)
approx      <- approx/sum(approx)

# Build plot data
pltData     <- data.frame( offset=offset, Exact=template, Gaussian=approx)
pltData     <- melt(pltData, id.vars="offset", variable_name="Method")

# Build plot
plt         <- ggplot(pltData, aes(x=offset, y=value, color=variable,
                                   fill=variable, shape=variable)) +
  geom_line() + geom_point() +
  xlab("Offset from nucleosome center") +
  ylab("Expected counts per nucleosome")
plt         <- plt + scale_colour_manual(values=c("red","black","blue"))
plt         <- plt + scale_fill_manual(values=c("red","black","blue"))

# Output plot
theme_set(theme_bw())

pdf(str_c(plotPathBase,".pdf"), width=6, height=6*0.62)
print(plt)
dev.off()

# Convert plot to eps
system(sprintf(convertCmdFmt, plotPathBase))
