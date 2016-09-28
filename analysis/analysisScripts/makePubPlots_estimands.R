library(ggplot2)
library(stringr)
library(data.table)

kSeed <- 20131109
kOutput <- "plots/pub/cluster_estimands.pdf"
kPlotSettings <- list(width=10, height=3.1, dpi=300)
kTheme <- theme_bw()
kLength <- 147
kR <- 0.5
kExamples <- data.frame(
  pos=rep(1:kLength, 4),
  # Reordering cases by concentration
  case=factor(rep(c("Single", "Uniform", "Symmetric", "Asymmetric"),
                  each=kLength),
              levels=c("Uniform", "Symmetric", "Asymmetric", "Single")),
  beta=c(
    # Centered spike
    c(rep(0, floor(kLength/2)), 1, rep(0, floor(kLength/2))),
    # Uniform
    rep(1, kLength),
    # Symmetric
    c(rep(0, 29), kR^4, rep(0, 10), kR^3, rep(0, 10), kR^2,
      rep(0, 10), kR, rep(0, 10), 1, rep(0, 10), kR,
      rep(0, 10), kR^2, rep(0, 10), kR^3, rep(0, 10), kR^4, rep(0, 29)) +
      rexp(kLength) * 0.005,
    # Skewed
    c(rep(0, 9), kR, rep(0, 15), kR^4, rep(0, 15), kR^2,
      rep(0, 15), kR^3, rep(0, 15), 1, rep(0, 10), kR^2,
      rep(0, 10), kR, rep(0, 10), kR^3, rep(0, 10), kR^4, rep(0, 29)) + 
      rexp(kLength) * 0.005
  )
)

ClusterEstimands <- function(beta, q=0.9) {
  p <- beta / sum(beta)
  n <- length(beta)
  k <- seq(n)
  # Center is the average position within the cluster
  center <- sum(p * k)
  # Localization is the mean average deviation of positions from the center
  localization <- 1 - 4 * mean(p * abs(k - center))
  # Structure is the normalized entropy of positions in the cluster
  structure <- 1 + sum(p * log(p + .Machine$double.eps)) / log(n)
  # Sparsity is the normalized quantile of positions in the cluster
  nq <- min(which(cumsum(sort(p, decreasing=TRUE)) >= q))
  sparsity <- 1 - (nq - 1) / q / n
  return(list(center=center,
              localization=localization,
              structure=structure,
              sparsity=sparsity))
}

EstimandsByCase <- function(.df, key="case") {
  data.table(.df, key=key)[, ClusterEstimands(beta), keyby=key]
}

BuildDiagram <- function(obs, ...) {
  estimands <- EstimandsByCase(obs, ...)
  estimands$label <- sprintf(paste0("paste(L == %0.1f, '; ', ",
                                    "S == %0.1f, '; ', ",
                                    "R == %0.1f)"),
                             estimands$localization,
                             estimands$structure,
                             estimands$sparsity)

  plt <- ggplot(data=obs, mapping=aes(x=pos, y=beta)) +
    geom_step() +
    facet_grid(. ~ case) +
    ylab(expression(beta[k])) + xlab("Position within cluster") +
    xlim(c(0, kLength + 1)) +
    geom_text(data=estimands, aes(label=label), x=kLength / 2, y=1.1,
              size=5, family="serif", parse=TRUE) +
    ylim(c(0, 1.1)) +
    kTheme
  return(plt)
}

SaveDiagram <- function(plt, file, ...) {
  ggsave(file, plot=plt, ...)
}

Main <- function() {
  set.seed(20131109)
  plt <- BuildDiagram(kExamples)
  do.call(SaveDiagram, c(list(plt, kOutput), kPlotSettings))
}

if (sys.nframe() == 0) {
  Main()
}