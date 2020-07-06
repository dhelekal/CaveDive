#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double coalescent_loglh(NumericVector time_intervals,
                                     NumericVector lineage_count,
                                     double pop_size
                                      ) {

  int n = time_intervals.length();
  auto theta = pop_size;
  double log_lh = 0;

  for(int i=0; i<n; ++i) {
      log_lh += -time_intervals[i]*lineage_count[i]*(lineage_count[i]-1)/(2*theta);
  }
  log_lh += -((n+1)/2-1)*std::log(theta);
  return(log_lh);
}

// [[Rcpp::export]]
double exponential_coalescent_loglh(NumericVector sampling_times,
                                     NumericVector coalescent_times,
                                     double pop_size,
                                     double lambda
                                     ) {
  auto s_idx = 0;
  auto c_idx = 0;

  double t = 0;
  double t_max = sampling_times[s_idx];
  int k = 1;

  double log_lh = 0;

  while (c_idx < coalescent_times.length()) {
    double delta_t;
    if (s_idx+1 < sampling_times.length() && sampling_times[s_idx+1] < coalescent_times[c_idx]){
        ++s_idx;
        delta_t = t_max - t - sampling_times[s_idx];
        log_lh += lambda * (t+delta_t) - k*(k-1)/(2*pop_size*lambda)*(std::exp(lambda*(t+delta_t))-std::exp(lambda*t));
        ++k;
        t += delta_t;
    } else {
        delta_t = t_max - t - coalescent_times[c_idx];
        ++c_idx;
        log_lh += lambda * (t+delta_t) - k*(k-1)/(2*pop_size*lambda)*(std::exp(lambda*(t+delta_t))-std::exp(lambda*t));
        --k;
        t += delta_t;
    }
  }

  log_lh += -(sampling_times.length()-1)*std::log(pop_size);
  return(log_lh);
}
