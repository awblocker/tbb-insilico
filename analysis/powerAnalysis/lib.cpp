#include <Rcpp.h>
using namespace Rcpp;

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

