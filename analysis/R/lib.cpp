#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector LabelClusters(const IntegerVector centers, const int n,
                            const int w) {
  /*
  Label positions with cluster center for cluster of width 2w + 1.
  Note that centers are expected to use 1-based indexing.
  */
  int n_clusters = centers.size();
  IntegerVector cluster_centers(n);
  int cluster=0, i=0, start_cluster, end_cluster;
  
  for (cluster=0; cluster < n_clusters; cluster++) {
    start_cluster = std::max(0, centers[cluster] - w - 1);
    end_cluster = std::min(n, centers[cluster] + w);
    for (i = start_cluster; i < end_cluster; i++) {
     cluster_centers[i] = centers[cluster]; 
    }
  }
  
  return cluster_centers;
}

// [[Rcpp::export]]
DataFrame ComputeClusterIndices(const IntegerVector centers,
                                const NumericVector x,
                                const int w, const double adj,
                                const double q) {
  /*
  Compute cluster-level summaries for clusters of width 2w + 1 with centers
  from centers vector based on statistics in x.
  
  Currently computes structure (entropy-based), localization (MAD-based),
  and sparsity indices.
  */
  int n_clusters = centers.size(), n_positions = x.size();
  NumericVector structure(n_clusters), localization(n_clusters),
                sparsity(n_clusters);
  NumericVector p(2 * w + 1), p_sorted(2 * w + 1), p_cumulative(2 * w + 1);
  int cluster=0, i=0, start_cluster, end_cluster, length_cluster;
  double sum_cluster, mean_abs_dev, entropy, center, p_regularized;
  
  for (cluster=0; cluster < n_clusters; cluster++) {
    // Get limits of cluster
    start_cluster = std::max(0, centers[cluster] - w - 1);
    end_cluster = std::min(n_positions, centers[cluster] + w);
    length_cluster = end_cluster - start_cluster;
    
    // Compute sum of entries within cluster
    sum_cluster = 0;
    for (i=0; i < length_cluster; i++) {
      sum_cluster += x[start_cluster + i];
    }
    
    // Compute entropy and center of cluster
    center = 0;
    entropy = 0;
    for (i = 0; i < length_cluster; i++) {
      p[i] = x[start_cluster + i] / sum_cluster;
      p_regularized = (x[start_cluster + i] + adj) / (sum_cluster + adj *
              length_cluster);
      entropy += -p_regularized * log(p_regularized);
      center += i * p[i];
    }

    // Compute mean absolute deviation
    mean_abs_dev = 0;
    for (i = 0; i < length_cluster; i++) {
      mean_abs_dev += p[i] * abs(i  - center);
    }
    
    // Create sorted version of p
    std::copy(p.begin(), p.begin() + length_cluster, p_sorted.begin());
    std::sort(p_sorted.begin(), p_sorted.begin() + length_cluster,
        std::greater<double>());

    // Compute sparsity from sorted version
    std::partial_sum(p_sorted.begin(), p_sorted.begin() + length_cluster,
        p_cumulative.begin());    
    for (i = 0; i < length_cluster; i++) {
      if (p_cumulative[i] > q) {
        break;
      }
    }
    sparsity[cluster] = 1. - (double) i / (length_cluster * q);

    
    // Compute and store normalized versions
    structure[cluster] = 1. - entropy / log(length_cluster);
    localization[cluster] = 1. - 4 * mean_abs_dev / length_cluster;
  }
  
  // Build dataframe of return values
  DataFrame out;
  out = DataFrame::create(Named("center") = centers,
                          Named("structure") = structure,
                          Named("localization") = localization,
                          Named("sparsity") = sparsity);
  return out;
}
// [[Rcpp::export]]
int GreedySearch(IntegerVector candidates, int min_distance,
                 IntegerVector peaks) {
  int n_candidates = candidates.size();
  int n_peaks = 1, candidate_distance = 0;
  
  // Start with first candidate
  peaks[0] = candidates[0];
  
  if (n_candidates > 1) {
    int k;
    bool accept;
    
    for (int i=1; i < n_candidates; i++) {
      // Compute minimum distance of next candidate to previous peaks
      accept = true;
      for (k=0; k < n_peaks; k++){
        if (abs(candidates[i] - peaks[k]) < min_distance) {
          accept = false;
          break;
        }
      }
      
      // Skip to next candidate if too close to previously called peak
      if (accept) {
        // Otherwise, add candidate and its value to the list
        n_peaks++;
        peaks[n_peaks - 1] = candidates[i];
      }
    }
  }
  
  return n_peaks;
}
