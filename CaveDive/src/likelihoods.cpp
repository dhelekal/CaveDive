#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double coalescent_likelihood(NumericVector time_intervals,
                                     NumericVector lineage_count,
                                     double pop_size
                                      ) {

  int n = time_intervals.length();
  double theta = pop_size;
  double log_lh = 0;

  for(int i=0; i<n; ++i) {
      log_lh += -time_intervals[i]*lineage_count[i]*(lineage_count[i]-1)/(2*theta);
  }
  log_lh += -((n+1)/2-1)*std::log(theta);
  return(log_lh);
}
