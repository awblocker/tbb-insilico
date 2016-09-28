# Load libraries
library(ggplot2)
library(grid)
library(gridExtra)
library(RcppCNPy)

source('powerAnalysis/lib.R')

# Constants

kDesignPath <- 'powerAnalysis/data/simChrom_design.txt'
kRepGenePath <- 'powerAnalysis/data/representative_gene.txt'
kSimYPath <- 'powerAnalysis/data/simChrom_y.txt'
kSimBetaPath <- 'powerAnalysis/data/simChrom_b.npy'
kGeneSpec <- c(coverage_quantile=0.55, alt_pos_spacing=25,
               alt_pos_magnitude=0.5, rep=1)
kGeneLength <- 3501
kPromoterLength <- 1000

kGeneStructurePlotPath <- 'powerAnalysis/plot_repGene.pdf'
kGeneReadsPlotPath <- 'powerAnalysis/plot_repGene_y.pdf'
kGeneCoefPlotPath <- 'powerAnalysis/plot_repGene_b.pdf'
kGeneCombinedPlotPath <- 'powerAnalysis/plot_repGene_combined.pdf'

# Load design data
design <- read.csv(kDesignPath)
design$coverage.bin <- as.integer(as.factor(design$coverage))

# Load data
repGene <- read.csv(kRepGenePath)
y.df <- ReadRaggedArray(kSimYPath, integer(0),
                        sep=',', value.var='y', row.var='rep')
b.matrix <- npyLoad(kSimBetaPath, dotranspose=TRUE)

# Extract gene of interest
gene.id <- with(
  design, which(coverage_quantile == kGeneSpec['coverage_quantile'] &
                  alt_pos_spacing == kGeneSpec['alt_pos_spacing'] & 
                  alt_pos_magnitude == kGeneSpec['alt_pos_magnitude']))
y.gene <- y.df$y[y.df$rep == kGeneSpec['rep']][
  (kGeneLength * (gene.id - 1) + 1) : (kGeneLength * gene.id)]
b.gene <- b.matrix[kGeneSpec['rep'],
                   (kGeneLength * (gene.id - 1) + 1) : (kGeneLength * gene.id)]

gene.df <- data.frame(pos=seq(kGeneLength) - kPromoterLength - 1,
                      y=y.gene, b=b.gene)


# Build plots

# Representative gene organization
plot.structure <- qplot(x=center, y=relative_occupancy, xend=center, yend=0,
                        data=repGene, geom='segment',
                        xlab='Position relative to TSS',
                        ylab='Relative occupancy')

# Reads
plot.y <- qplot(x=pos, y=y, xend=pos, yend=0, data=gene.df, geom='segment',
                xlab='Position relative to TSS', ylab=expression(y[k]))

# Coefficients
plot.b <- qplot(x=pos, y=b, xend=pos, yend=0, data=gene.df, geom='segment',
                xlab='Position relative to TSS', ylab=expression(beta[k]))

# Output plots
theme_set(theme_bw())

pdf(kGeneStructurePlotPath, width=5, height=3.1)
print(plot.structure)
dev.off()

pdf(kGeneReadsPlotPath, width=5, height=3.1)
print(plot.y)
dev.off()

pdf(kGeneCoefPlotPath, width=5, height=3.1)
print(plot.b)
dev.off()

pdf(kGeneCombinedPlotPath, width=10, height=4)
grid.arrange(plot.y + theme(axis.title.x=element_blank()),
             plot.b + theme(axis.title.x=element_blank()),
             nrow=2, ncol=1,
             sub='Position relative to TSS')
dev.off()
