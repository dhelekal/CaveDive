#include <Rcpp.h>
using namespace Rcpp;

double half_log_int(double r, double N, double t0, double t, double delta_t){

  double out = std::numeric_limits<double>::infinity();
  if (t+t0 < 0 && t+delta_t+t0 < 0) {
    out =  ((-1/(r*N))*(std::log(std::exp(-r*(t0+t+delta_t))-1)-std::log(std::exp(-r*(t0+t))-1)));
  } 

  return out;
}

double half_log_rate(double r, double N, double t0, double t){
  return 1/(N*(1-exp(r*(t+t0))));
}

double const_int(double N, double t0, double delta_t){
  return delta_t/N;
}

double const_rate(double N, double t){
  return 1/N;
}

double exp_int(double lambda, double N, double t, double delta_t) {
    return (1/(N*lambda))*(std::exp(lambda*(t+delta_t))-std::exp(lambda*t));
}

double exp_rate(double lambda, double N, double t) {
    return (1/N) * std::exp(lambda*t);
}

//'Compute likelihood of a realisation of a homogeneous coalescent.
//' 
//' @param time_intervals waiting times.
//' @param lineage_count number of lineages for each waiting time.
//' @param pop_size effective population size \code{Neg}.
//' @return the log likelihood for the given parameters.
//' @export
// [[Rcpp::export]]
double coalescent_loglh(NumericVector sampling_times,
                        NumericVector coalescent_times,
                        double pop_size) {
    auto s_idx = 0;
    auto c_idx = 0;

    double t = 0;
    double t_max = sampling_times[s_idx];
    int k = 1;

    double log_lh = 0;

    while (c_idx < coalescent_times.length()) {
        double delta_t;
        if (s_idx+1 < sampling_times.length() && coalescent_times[c_idx] < sampling_times[s_idx+1]){
            //Sampling
            ++s_idx;
            delta_t = t_max - t - sampling_times[s_idx];
            log_lh += (-k*(k-1)/(2))*const_int(pop_size, t, delta_t);
            ++k;
            t += delta_t;
        } else {
            //Coalescent
            delta_t = t_max - t - coalescent_times[c_idx];
            ++c_idx;
            log_lh += std::log(const_rate(pop_size,t+delta_t)) + (-k*(k-1)/(2))*const_int(pop_size, t, delta_t);
            --k;
            t += delta_t;
        }
    }
    return(log_lh);
}

//'Compute likelihood of a realisation of an exponential growth coalescent.
//' 
//' @param sampling_times times of leaves.
//' @param coalescent_times times of coalescent events.
//' @param lambda growth rate.
//' @param pop_size final effective pop_size at time of most recent sample.
//' @return the log likelihood for the given parameters.
//' @export
// [[Rcpp::export]]
double exponential_coalescent_loglh(NumericVector sampling_times,
                                     NumericVector coalescent_times,
                                     const double lambda,
                                     const double pop_size) {

  auto s_idx = 0;
  auto c_idx = 0;

  double t = 0;
  double t_max = sampling_times[s_idx];
  int k = 1;

  double log_lh = 0;

  while (c_idx < coalescent_times.length()) {
    double delta_t;
    if (s_idx+1 < sampling_times.length() && coalescent_times[c_idx] < sampling_times[s_idx+1]){
        ++s_idx;
        delta_t = t_max - t - sampling_times[s_idx];
        log_lh += (-k*(k-1)/(2))*(exp_int(lambda, pop_size, t, delta_t));
        ++k;
        t += delta_t;
    } else {
        delta_t = t_max - t - coalescent_times[c_idx];
        ++c_idx;
        log_lh += std::log(exp_rate(lambda, pop_size, t+delta_t)) + (-k*(k-1)/(2))*(exp_int(lambda, pop_size, t, delta_t));
        --k;
        t += delta_t;
    }
  }

  return(log_lh);
}

//'Compute likelihood of a realisation of an half/logistic growth coalescent.
//' 
//' @param sampling_times times of leaves.
//' @param coalescent_times times of coalescent events.
//' @param t0 time of divergence.
//' @param growth rate.
//' @param N asymptotic population size.
//' @return the log likelihood for the given parameters.
//' @export
// [[Rcpp::export]]
double logexp_coalescent_loglh(NumericVector sampling_times,
                                     NumericVector coalescent_times,
                                     const double t0,
                                     const double r,
                                     const double N
                                     ) {
  auto s_idx = 0;
  auto c_idx = 0;

  double t = 0;
  double t_max = sampling_times[s_idx];
  int k = 1;

  double log_lh = 0;

  while (c_idx < coalescent_times.length()) {
    double delta_t;
    if (s_idx+1 < sampling_times.length() && coalescent_times[c_idx] < sampling_times[s_idx+1]){
        //Sampling
        ++s_idx;
        delta_t = t_max - t - sampling_times[s_idx];
        log_lh += (-k*(k-1)/(2))*half_log_int(r, N, t0, t, delta_t);
        ++k;
        t += delta_t;
    } else {
        //Coalescent
        delta_t = t_max - t - coalescent_times[c_idx];
        ++c_idx;
        log_lh += std::log(half_log_rate(r, N, t0, t+delta_t)) + (-k*(k-1)/(2))*half_log_int(r, N, t0, t, delta_t);
        --k;
        t += delta_t;
    }
  }
  return(log_lh);
}

